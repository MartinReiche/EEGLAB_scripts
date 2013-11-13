% Find number of given triggers with given ERP window in EEG file and extract
% epochs using these triggers
%
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

function eventCount = segmentation(EEG,trig,analysis,paths,eventCount,iFile,iSubj,condOrder)
    
% count all events in current EEG file
    disp(':: Get events for current file');
    allEvents = [];
    % Switch rejection Mode
    switch analysis.rejmode
        % Delta rejection only
      case {1,3,4}
        
        % get vector of all events
        for iEvent = 1:size(EEG.event,2)
            allEvents = [allEvents {EEG.event(iEvent).type}];
        end

        % Get numbers of all Events
        for iTrig = 1:size(trig.triggers,1)
            eventCount(iTrig,iFile) = sum(ismember(allEvents,trig.triggers{iTrig,1}));
        end
        
        % start epoching available triggers of current block
        for iTrig = 1:size(trig.triggers,1)
            if eventCount(iTrig,iFile) > 0 
                % if there are available events for the current trigger type
                disp([':: Extracting segments for Trigger ' num2str(trig.triggers{iTrig,1})]);
                disp(' ');
                
                currTrig = trig.triggers{iTrig,1};
                EEG2 = pop_epoch(EEG,{currTrig},[0.001*analysis.erpWin(1) ...
                                    0.001*analysis.erpWin(2)]);
                % if rejection mode is sorted averaging
                if analysis.rejmode == 4
                    % perform baseline correction
                    EEG2 = pop_rmbase(EEG2,analysis.baseWin);
                end
                % save trigger set for current block
                EEG2 = pop_saveset(EEG2,[paths.resFileSubSpec num2str(iSubj,'%0.2d')...
                                    paths.resFileBlockSpec num2str(iFile, '%0.2d') ...
                                    paths.resFileTrigSpec num2str(trig.triggers{iTrig,1},'%0.2d') ...
                                    '.set'],paths.resDir);
                disp(' ');
            end
        end
        % delta rejection plus EMC
      case 2
        clear EEG;
        EEG = pop_loadset([paths.resDir...
                           paths.resFileSubSpec num2str(iSubj, '%0.2d')...
                           'all.set']);        

        allEvents = [];
        
        % Get vector of all events
        for iEvent = 1:size(EEG.event,2)
            allEvents = [allEvents {EEG.event(iEvent).type}];
        end                
        
        % Get numbers of all Events
        eventCount = zeros(size(trig.triggers,1),1);   

        for iSubTrig = 1:size(trig.triggers,1)
            eventCount(iSubTrig) = sum(ismember(allEvents,trig.triggers{iSubTrig,1}));
        end

        % start epoching availabele trig.triggers of current block
        for iTrig = 1:size(trig.triggers,1)
            if eventCount(iTrig) > 0
                %if there are available events for the current trigger type
                disp([':: Extracting segments for Trigger ' num2str(trig.triggers{iTrig,1})]);
                disp(' ');
                % convert current triggers to cell array of strings
                currTrig = trig.triggers{iTrig,1};
                EEG2 = pop_epoch(EEG,{currTrig},[0.001*analysis.erpWin(1) ...
                                    0.001*analysis.erpWin(2)]);
                % save EEG
                EEG2 = pop_saveset(EEG2,[paths.resFileSubSpec num2str(iSubj,'%0.2d')...
                                    paths.resFileTrigSpec num2str(trig.triggers{iTrig,1},'%0.2d')...
                                    'all.set'],paths.resDir);
                clear EEG2;
                disp(' ');
            else
                disp([':: No available triggers found for event ' ...
                     num2str(trig.triggers{iTrig,1}) ', skipping.']);
            end
        end
                
        % save numbers of events in result dir of current subject
        disp(':: Saving initial number of events');
        save([paths.resDir paths.resFileSubSpec num2str(iSubj, '%0.2d') 'nTrials.mat'],'eventCount');
        % delete old epoching all data file        
        disp(':: Deleting old all data file');
        delete([paths.resDir paths.resFileSubSpec num2str(iSubj, '%0.2d') 'all.*']);
    end
