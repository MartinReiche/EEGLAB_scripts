function EEG = checkFileBasic(EEG,nSubj,iFile,taskType,analysis,trig,paths,condOrder)
% Check the currently loaded raw data file for sampling rate, triggercodes
% and trigger latencies
%
% Copyright (c) 2013 Martin Reiche, University of Leipzig
% Author: Martin Reiche, reiche.stud@gmail.com

% Random Changes
EEG = subjChanges(EEG,nSubj,iFile,taskType,analysis,trig,paths,condOrder);

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
if strcmpi(analysis.rawFormat,'biosemi')
    % convert all triggers to strings BioSemi Raw Data
    disp(' ');
    disp(':: Converting all triggers to strings');
    for iTrig = 1:size(EEG.event,2)
        EEG.event(iTrig).type = num2str(EEG.event(iTrig).type);
        EEG.urevent(iTrig).type = num2str(EEG.urevent(iTrig).type);
    end

end

% remove unintended events
for iDel = 1:size(trig.delete,2)
    iTrig = 1;
    while iTrig <= size(EEG.event,2)
        if strcmpi(EEG.event(iTrig).type,trig.delete{iDel})
            disp([':: Removing ' trig.delete{iDel} ' event']);
            EEG.event(iTrig) = [];
        else
            iTrig = iTrig + 1;
        end
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

%% Check for trigger types
if taskType == 1
    for iTrig = 1:numel(eegEvent)
        if ~any(strcmp(trig.checkTrig{condOrder(iFile)},eegEvent(iTrig))); 
            disp([':: Detected wrong trigger (' eegEvent{iTrig} ') at position: ' num2str(iTrig) ' for Subject '...
                  num2str(nSubj) ' Block ' num2str(iFile)]);
            input(':: Press Ret to continue, C-c to abort.');
        end
    end
end
% % check number of triggers
% if numel(EEG.event) ~= trig.num
%     disp(['Trigger number mismatch in File: ' num2str(iFile) ' of subject ' num2str(nSubj)]);
%     disp([':: expected number of triggers: ' num2str(trig.num)]);
%     disp([':: detected number of triggers: ' num2str(numel(EEG.event))]);
%     input(':: Press Ret to continue, C-c to abort.');
% else
%     disp([':: Checked Trigger for Subject ' num2str(nSubj) ' Block ' num2str(iFile) '. Success.']);
% end

%% Check trigger timing
if taskType == 1
    % get intended trigger onset time
    eegLat = diff(eegTime);
    % convert to milliseconds
    eegLat = abs((eegLat .* checkFactor) - (trig.soa * 1000));

    if sum(eegLat > 2) >= trig.rejThresh.stim; 
        % mark missed triggers
        for iTrig = 2:numel(eegLat)
            if eegLat(iTrig) > 2
                EEG.event(iTrig).type = num2str(trig.missTrig);
            end
        end
    end

    % save number of rejected triggers
    eventExcp = sum(eegLat > 2);

    if sum(eventExcp) > trig.rejThresh.note
        disp(' ');
        disp([':: More than ' num2str(trig.rejThresh.note) ' Events have been rejected for Subject '...
              num2str(nSubj) ' Block ' num2str(iFile)]);
        disp([':: Number of rejected events: ' num2str(eventExcp)]);
        if trig.rejThresh.pause && ~strcmp(paths.cluster,'remote')
            input(':: Press Ret to continue, C-c to abort.');
        end
    end
end