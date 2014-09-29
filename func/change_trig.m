% Changes triggers at start of block to exclude them from the analysis,
% change omission triggers to hit vs. miss omissions and exclude events in
% predefined ranges around triggers and response trials

% Copyright (c) 2013 Martin Reiche, Carl-von-Ossietzky-University Oldenburg
% Author: Martin Reiche, martin.reiche@uni-oldnburg.de

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


function EEG = change_trig(EEG,analysis,trig,iFile,condOrder,taskType)
% first trigger is already recoded, also recode the second trigger of each
% block that it won't be analyzed

disp(':: Retrigger events at block start');
for iTrig = 1:numel(analysis.eventExcl);
    % if trigger is not already excluded
    if ~any(ismember(EEG.event(analysis.eventExcl(iTrig)).type,trig.exclBefore))
        EEG.event(analysis.eventExcl(iTrig)).type = analysis.exclTrig;
    end
end

% define trigger range of omissions

%% Exclusion of trials in predefined ranges around omissions and responses
for iTrig = 1:size(EEG.event,2)    
    % go through all events in current block and check for response triggers
    if strcmp(EEG.event(iTrig).type,analysis.respTrig)
        %% IF RESPONSE TRIGGER IS ENCOUNTERED
        iPreTrig = iTrig;
        iPostTrig = iTrig;
        % Correct all trials BEFORE response
        searchRange = 1;
        while iPreTrig > 1 && searchRange
            % count 1 back from last response
            iPreTrig = iPreTrig - 1;
            % get latency range from last response
            preRange = ((EEG.event(iPreTrig).latency / analysis.sampRate) - (EEG.event(iTrig).latency / analysis.sampRate)) * 1000;
            if preRange >= analysis.respEx(1)
                % if the current event is not an omission or a response
                if ~ismember(EEG.event(iPreTrig).type,analysis.omissionRange) |  ~strcmp(EEG.event(iPreTrig).type,analysis.respTrig)
                    EEG.event(iPreTrig).type = analysis.exclTrig;
                end
            else
                searchRange = 0; 
            end
        end
        % Correct all trials AFTER response
        searchRange = 1;
        while iPostTrig < size(EEG.event,2) && searchRange
            % add 1 from last response
            iPostTrig = iPostTrig + 1;
            % get latency range from last response
            postRange = ((EEG.event(iPostTrig).latency / analysis.sampRate) - (EEG.event(iTrig).latency / analysis.sampRate)) * 1000;
            if postRange <= analysis.respEx(2)
                % if the current event is not an omission or a response
                if ~ismember(EEG.event(iPostTrig).type,analysis.omissionRange) | ~strcmp(EEG.event(iPostTrig).type,analysis.respTrig)
                    EEG.event(iPostTrig).type = analysis.exclTrig;
                end
            else
                searchRange = 0; 
            end
        end
    end
    if ismember(EEG.event(iTrig).type,analysis.omissionRange)
        %% IF OMISSION TRIGGER IS ENCOUNTERED
        iPreTrig = iTrig;
        iPostTrig = iTrig;
        % Correct all trials BEFORE response
        searchRange = 1;
        while iPreTrig > 1 && searchRange
            % count 1 back from last response
            iPreTrig = iPreTrig - 1;
            % get latency range from last response
            preRange = ((EEG.event(iPreTrig).latency / analysis.sampRate) - (EEG.event(iTrig).latency / analysis.sampRate)) * 1000;
            if preRange >= analysis.omEx(1)
                % if the current event is not an omission or a response
                if ~ismember(EEG.event(iPreTrig).type,analysis.omissionRange) || ~strcmp(EEG.event(iPreTrig).type,analysis.respTrig)
                    EEG.event(iPreTrig).type = analysis.exclTrig;
                end
            else
                searchRange = 0; 
            end
        end
        % Correct all trials AFTER response
        searchRange = 1;
        while iPostTrig < size(EEG.event,2) && searchRange
            % add 1 from last response
            iPostTrig = iPostTrig + 1;
            % get latency range from last response
            postRange = ((EEG.event(iPostTrig).latency / analysis.sampRate) - (EEG.event(iTrig).latency / analysis.sampRate)) * 1000;
            if postRange <= analysis.omEx(2)
                % if the current event is not an omission or a response
                if ~ismember(EEG.event(iPostTrig).type,analysis.omissionRange) || ~strcmp(EEG.event(iPostTrig).type,analysis.respTrig)
                    EEG.event(iPostTrig).type = analysis.exclTrig;
                end
            else
                searchRange = 0; 
            end
        end
    end
end

%% Retrigger Omission (Hit vs Miss)
if taskType == 2

    % in position 1, 70 dB SPL:
    % hit = 123 followed by 1
    % miss = 123 followed by 2
    % false alarm = 113 followed by 2
    % correct rejection = 113 followed by 1

    % in position 2, 70 dB SPL:
    % hit = 124 followed by 1
    % miss = 124 followed by 2
    % false alarm = 114 followed by 2
    % correct rejection = 114 followed by 1

    for iEvent = 1:numel(EEG.event)
        % go through all events
        if ismember(EEG.event(iEvent).type,{'123', '113', '124', '114'})
            foundResponse = 0;
            iPostResponse = iEvent;
            noiseTrig = EEG.event(iEvent).type;
            
            while ~foundResponse
                % add one to increment variable
                iPostResponse = iPostResponse + 1;
                
                if ismember(EEG.event(iPostResponse).type,{'1', '2', '3', '4'})
                    foundResponse = 1;
                    responseTrig = EEG.event(iPostResponse).type;
                    % Retrigger noise stimuli according to response
                    if responseTrig == '1' 
                        if ismember(noiseTrig,{'123', '124'})
                            % HIT
                            EEG.event(iEvent).type = num2str(str2num(EEG.event(iEvent).type) + 1000);
                        elseif ismember(noiseTrig,{'113', '114'})
                            % CORRECT REJECTION
                            EEG.event(iEvent).type = num2str(str2num(EEG.event(iEvent).type) + 2000);
                        end
                    elseif responseTrig == '2'
                        if ismember(noiseTrig,{'123', '124'})
                            % MISS
                            EEG.event(iEvent).type = num2str(str2num(EEG.event(iEvent).type) + 3000);
                        elseif ismember(noiseTrig,{'113', '114'})
                            % FALSE ALARM
                            EEG.event(iEvent).type = num2str(str2num(EEG.event(iEvent).type) + 4000);
                        end
                    end
                end
            end
        end
    end
end

% change trigger latencies if intended
if ~isempty(trig.changeLat)
   % calculating latency change in sampling points
   latShift = round((analysis.sampRate/1000)*trig.changeLat);
   disp([':: Correcting trigger latencies. Shift latencies by ' num2str(trig.changeLat) ' ms (' num2str(latShift) ' sampling points)']);
   % shifting latency for all triggers
   for iTrig = 1:size(EEG.event,2)
       EEG.event(iTrig).latency =  EEG.event(iTrig).latency + latShift;
   end
end

% change the triggers as indiated by func/retrigConf.m
if analysis.changeTrig 
    reTrig = triggerlabels('retrig',trig,taskType,condOrder(iFile));
    % retrigger events
    disp(':: Retriggering');
    % start retriggering according to specified triggers in retrigConf.m
    % leave first and last trigger of file untouched
    for iEvent = 2:size(EEG.event,2)-1
        for iReTrig = 1:size(reTrig,1)
            if ismember(EEG.event(iEvent).type,reTrig{iReTrig,1})
                EEG.event(iEvent).type = reTrig{iReTrig,2};
            end
        end
    end
end
