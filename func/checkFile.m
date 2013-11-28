function [EEG, eventExcp] = checkFile(EEG,nSubj,iFile,stimPars,taskType,analysis,trig)
% Check the currently loaded raw data file for sampling rate, triggercodes
% and trigger latencies
%
% Copyright (c) 2013 Martin Reiche, University of Leipzig
% Author: Martin Reiche, reiche.stud@gmail.com

%% Parameters
checkFactor = 1000 / analysis.sampRate;




%% check sampling rate
if EEG.srate ~= analysis.sampRate
    error([...
        ':: Wrong Sampling Rate detected for Subject '...
        num2str(nSubj) ' in Block ' num2str(iFile) '. ('...
        num2str(EEG.srate) ' Hz)']);
end

%% Check Trigger
switch analysis.rawFormat
  case 'biosemi'
    % convert all triggers to strings BioSemi Raw Data
    disp(':: Converting all triggers to strings');
    for iTrig = 1:size(EEG.event,2)
        EEG.event(iTrig).type = num2str(EEG.event(iTrig).type);
        EEG.urevent(iTrig).type = num2str(EEG.urevent(iTrig).type);
    end
    % convert triggers to strings from stimulation file
    trigger = cell(size(stimPars.pars,2),1);
    for iTrig = 1:size(stimPars.pars,1)
       trigger{iTrig} = num2str(stimPars.pars(iTrig,5));
    end
  case 'brainvision'
    % convert trigger numbers to strings
    trigger = cellstr(num2str(stimPars.pars(:,5)));
    % add 'S' before all trigger strings from behavioral data
    for iTrig = 1:size(trigger,1)
        trigger{iTrig} = ['S' trigger{iTrig}];
    end
end

% get trigger timing
triggerTime = stimPars.pars(:,6);

% remove boundary events
iTrig = 1;
while iTrig <= size(EEG.event,2)
    if strcmpi(EEG.event(iTrig).type,'boundary')
       disp(':: Removing boundary event');
       EEG.event(iTrig) = [];
    else
        iTrig = iTrig + 1;
    end
end

% get EEG events and latency
eegEvent = cell(size(EEG.event,2),1);
% for event latency storage
eegTime = zeros(size(EEG.event,2),1);

for iTrig = 1:size(EEG.event,2)
    % disp([':: Loop: ' num2str(iTrig)])
    % disp([':: eegEvent: ' num2str(eegEvent(iTrig))])
    % disp([':: EEG.event: ' num2str(EEG.event(iTrig).type)])
    eegEvent{iTrig} = EEG.event(iTrig).type;
    eegTime(iTrig) = EEG.event(iTrig).latency;
end

% check for start and end trigger
if strcmpi(eegEvent{1},trig.startTrig)
    % if the start trigger is the first event, delete it
    disp(':: Removing start trigger')
    eegEvent{1} = [];
    eegTime(1) = [];
end
if strcmpi(eegEvent{end},trig.endTrig)
    % if the end trigger is the last event, delete it
    disp(':: Removing end trigger')
    eegEvent{end} = [];
    eegTime(end) = [];
end
% remove empty cells (where start end end trigger used to be)
eegEvent = eegEvent(~cellfun('isempty',eegEvent));

% ADDING RESPONSE TRIGGERS (RelAX, not configured for trigger strings)
% if task is active, add response triggers where intended
if taskType == 2
    trigger = trigger';
    triggerTime = triggerTime';
    % get all positions of responses
    iPos = find(stimPars.results(:,1) ~= 0)';
    % shift them for one position because we want a response trigger after
    % every event where a response took place
    iPos = iPos + 1;
    % define vector of responses
    respEvent = ones(1,numel(iPos)) .* trig.respTrig;
    respTime = stimPars.results(stimPars.results(:,1)~=0,1);
    % define new trigger array with new length: trigger + iPos
    newTrigger = zeros(1,length(trigger)+length(iPos));
    newTriggerTime = zeros(1,length(triggerTime)+length(iPos));
    % insert triggers at iPos
    newTrigger(iPos + (0:length(iPos)-1)) = respEvent;
    newTriggerTime(iPos + (0:length(iPos)-1)) = respTime;
    % insert old triggers
    newTrigger(~newTrigger) = trigger;
    newTriggerTime(~newTriggerTime) = triggerTime;
    trigger = newTrigger';
    triggerTime = newTriggerTime';
end

missCounter = 0;
% check whether there is a mismatch of triggers
if (size(trigger,1) ~= size(eegEvent,1))
    disp(' ')
    disp([':: Numbers of triggers do not match for Subject '...
          num2str(nSubj) ' Block ' num2str(iFile)]);
    disp([':: Expected: ' num2str(size(trigger,1)) ' triggers'...
          ', detected: ' num2str(size(eegEvent,1)) ' triggers.'...
         ]);
    disp(' ');
    %    input(':: Proceed to trigger correction? (Ret, C-c to abort)');
    disp(' ');
    
    for iTrial = 1:numel(trigger)
        % compare all triggers with all eegEvents stepwise
        if ~strcmpi(trigger{iTrial},eegEvent{iTrial}) % || ~changeFlag
            
            % check whether EEG has more or less than intended triggers
            switch 1
              case size(trigger,1) > size(eegEvent,1)
                % if Triggers are missing in the raw EEG file
                disp(':: Missing triggers in EEG file detected.');
                disp([':: Detected trigger mismatch at trial ' num2str(iTrial)]);
                disp([':: Expected trigger: ' trigger{iTrial}...
                      ', detected trigger: ' eegEvent{iTrial} '.']);
                disp(':: Trigger will be marked for rejection.')
                disp(' ');
                % insert miss trigger to check wheter the rest is in
                % place
                eegEvent = [eegEvent(1:(iTrial - 1)); trig.missTrig; eegEvent((iTrial):end)];
                % add 1 to missCounter for Triggers
                missCounter = missCounter + 1;
                % insert misstrigger at that position
                EEG2 = EEG;
                EEG2.event(iTrial).type = trig.missTrig;
                EEG2.event(iTrial).latency = ceil(((triggerTime(iTrial) - triggerTime(iTrial - 1))*1000) / checkFactor) + EEG2.event(iTrial - 1).latency;
                for iEvent = iTrial:(size(EEG.event,1)+1);
                    EEG2.event(iEvent).type = EEG.event(iEvent+1).type;
                    EEG2.event(iEvent).latency = EEG.event(iEvent+1).latency;
                end
                EEG = EEG2;
                trigger(iTrial) = trig.missTrig;
                eegTime = [eegTime(1:(iTrial - 1)); EEG.event(iTrial).latency; eegTime((iTrial):end)];
                
              case size(trigger,1) < size(eegEvent,1)
                % if there are more trigger in the EEG Data, mark the
                % mismatching trigger in the EEG data for later rejection
                disp([':: Unintended trigger in EEG file detected at trial '...
                      num2str(iTrial) ', removing event (' eegEvent{iTrial}]);
                disp([':: Expected trigger: ' trigger{iTrial}...
                      ', detected trigger: ' eegEvent{iTrial} '.']);
                disp(' ');
                % mark unintended trigger for rejection
                % EEG.event(iTrig).type = trig.missTrig;
                % trigger = [trigger(1:(iTrial - 1)); trig.missTrig; trigger((iTrial):end)];
                EEG.event(iTrial) = [];
                eegEvent(iTrial) = [];
                eegTime(iTrial) = [];
                % remove empty cells (which had just been emptied)
                eegEvent = eegEvent(~cellfun('isempty',eegEvent));
            end
        end
    end
    %    input(':: Press Ret to proceed, C-c to abort.');
end
% if trigger numbers are equal between intended and detected
if all(strcmp(trigger,eegEvent))
    % everything as intended
    disp([':: Checked Trigger for Subject ' num2str(nSubj)...
          ' Block ' num2str(iFile) '. Success.']);
else
    % if trigger numbers match but trigger mismatch
    trig.missTrig = find(~strcmp(trigger,eegEvent));
    disp(' ');
    for iMsg = 1:size(trig.missTrig,1)
        disp([':: Trigger mismatch at trial: '...
              num2str(trig.missTrig(iMsg))...
              ' (expected: ' trigger{trig.missTrig(iMsg)}...
              ', got: ' eegEvent{trig.missTrig(iMsg)} ')'...
             ]);
    end
    disp(' ');
    error([':: Trigger mismatch detected for Subject ' num2str(nSubj)...
           ' Block ' num2str(iFile)]);
end

%% Check trigger timing
% get intended trigger onset time

eegLat = zeros((size(eegTime,1) - 1),1);
triggerLat = zeros((size(triggerTime,1) - 1),1);
for iTrial = 2:size(eegTime)
    % run through all triggers and get latency between two consecutive
    % triggers, if preceeding was response trigger, get the latency between
    % the stimulus trigger and the other stimulus trigger before
    if strcmp(trigger{iTrial-1},trig.respTrig)
        eegLat(iTrial) = eegTime(iTrial) - eegTime(iTrial - 2);
        triggerLat(iTrial) = triggerTime(iTrial) - triggerTime(iTrial - 2);
    else
        eegLat(iTrial) = eegTime(iTrial) - eegTime(iTrial - 1);
        triggerLat(iTrial) = triggerTime(iTrial) - triggerTime(iTrial - 1);
    end
end
% convert to milliseconds
eegLat = eegLat .* checkFactor;
eegLat = {trigger (triggerLat .*1000) eegLat};
eegLat = [eegLat abs(eegLat{:,3} - eegLat{:,2})];

% check if trigger latency exceeds predefined threshold
rejCounter = 0;

for iTrial = 1:size(eegLat,1) 
    if ((eegLat{iTrial,4} > trig.rejThresh.stim) & (eegLat{iTrial,1} ~= trig.respTrig)) | ((eegLat{iTrial,4} > trig.rejThresh.resp) & (eegLat{iTrial,1} == trig.respTrig)) 
        disp(':: Trigger latency exceeded accepted value, marked for rejection');
        disp([eegLat{iTrial,4} trig.rejThresh.stim]);
        EEG.event(iTrial).type = trig.missTrig;
        rejCounter = rejCounter + 1;
    end
end
% save missing and rejected events for output
eventExcp.rej = rejCounter;
eventExcp.miss = missCounter;
if sum(eventExcp.rej) > trig.rejThresh.note
    disp(' ');
    disp([':: More than ' num2str(trig.rejThresh.note) ' Events have been rejected for Subject '...
          num2str(nSubj) ' Block ' num2str(iFile)]);
    if analysis.rejTresh.pause
        input(':: Press Ret to continue, C-c to abort.');
    end
end
end
