% Perform artifact rejection, form average and return ERPs for given subject
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

function [subErp,subErpEqual,subTrialInd,corrTrials,eventCount,rejEpoch,trialNum] = eeg_rejection(iSubj,paths,trig,analysis)
% get current dir
subjDir = [paths.resDirAll paths.resSubFolderPrefix num2str(iSubj, '%0.2d') '/'];
% get list of  file names in current dir
subjFiles = dir([subjDir '*set'])
% load one (the 1st) of the .set files for the current subject
EEG = pop_loadset(subjFiles(1).name,subjDir);
% initialize erp matrix for curren subject
subErp = zeros(1,size(trig.triggers,1),size(EEG.data,2),size(EEG.data,1));
subErpEqual = zeros(1,size(trig.triggers,1),size(EEG.data,2),size(EEG.data,1));
% initialize cell arrays for vector of rejected epoch indices per trial
rejEpoch = cell(1,size(trig.triggers,1));
clear EEG;

% clear old reject files (just in case - nothing should be there)
if size(dir([subjDir '*good*']),1)
    % if some old epoch files already existed, delete them
    disp(' ');
    disp([':: Detected old corrected epoch files for subject ' ...
          num2str(iSubj,'%0.2d') ', deleting them']);
    delete([subjDir 'Subj' num2str(iSubj,'%0.2d') '*good*']);
    delete([subjDir 'Subj' num2str(iSubj,'%0.2d') '*bad*']);
end

% get number of maximal trials per trigger (eventCount, 3rd col)
loadData = load([subjDir 'Subj' num2str(iSubj, '%0.2d') ...
                'nTrials.mat'],'eventCount');
eventCount = loadData.eventCount;

% initialize folder for rejection information
if analysis.batchMode
    % rename ERP File
    paths.rejFolder = [paths.resDirAll 'JOB_' num2str(analysis.jobIndex) '_usable_Data_' paths.erpFileName(1:end-4) '/'];
    paths.erpFileName = ['JOB_' num2str(analysis.jobIndex) '_'  paths.erpFileName]; 
else
    paths.rejFolder = [paths.resDirAll 'usable_Data_' paths.erpFileName(1:end-4) '/'];
end

if ~exist(paths.rejFolder,'dir');
    disp(':: Creating new folder to store usable Data text files');
    mkdir(paths.rejFolder)
end
% open text file for report on rejected epoches
rejFile = fopen([paths.rejFolder 'Subj' num2str(iSubj, '%0.2d')...
                 'reject.txt'],'w');

%% find channel numbers to exclude for current Subject
sFound = 0;
s = 1;
exclVector = [];
while s <= size(analysis.excludeElecs,1) && ~sFound
    if strcmp(analysis.excludeElecs{s}(1,:),num2str(iSubj, '%0.2d'))
        sFound = 1;
        disp(' ');
        disp([':: Excluded electrodes for subject ' num2str(iSubj, '%0.2d') ':']);
        if strcmp(analysis.excludeElecs{s}(2,:),'none')
            disp('none');
            disp(' ');
        else
            for ne = 2:size(analysis.excludeElecs{s,:})
                disp(analysis.excludeElecs{s}(ne,:));
                eFound = 0;
                for iChannel = 1:size(EEG.data,1) % determine relevant channel numbers
                    if strcmp(EEG.chanlocs(iChannel).labels,analysis.excludeElecs{s}(ne,:))
                        eFound = 1;
                        exclVector = [exclVector iChannel];
                    end
                end
                if ~eFound
                    input(['electrode not found in data, go on anyway?']);
                end
            end
            disp([':: (channels ' num2str(exclVector) ')']);
            disp(' ');
        end
    else
        s = s+1;
    end
end
if ~sFound
    input(['subject ' num2str(iSubj, '%0.2d') ' not found in electrode exclusion list, go on including all their channels?']);
end % End of finding Channels to exclude

if ~isempty(exclVector)
    exclLabels = [];
    for iChannel = exclVector
        exclLabels = [exclLabels ' ' EEG.chanlocs(iChannel).labels];
    end
    fprintf(rejFile,'%s\n','----------------------------------------------------------------');
    fprintf(rejFile,'%s\n',['Rejection mode: ' analysis.rejLabel{analysis.rejmode+1,:} ', excluding channels' exclLabels]);
    fprintf(rejFile,'%s\n','----------------------------------------------------------------');
else
    fprintf(rejFile,'%s\n','----------------------------------------------------------------');
    fprintf(rejFile,'%s\n',['Rejection mode: ' analysis.rejLabel{analysis.rejmode+1,:} ', including all channels']);
    fprintf(rejFile,'%s\n','----------------------------------------------------------------');
end
clear EEG;

% initialize matrix for number of usable data per trigger
corrTrials = zeros(1,size(trig.triggers,1));

%% Start Rejection for all files
for iTrig = 1:size(trig.triggers,1)

    % get current file
    filename = [subjDir paths.resFileSubSpec num2str(iSubj, '%0.2d')...
                paths.resFileTrigSpec num2str(trig.triggers{iTrig,1},'%0.2d') 'all.set'];
    % if the current epoch file exists
    if exist(filename,'file')
        EEG = pop_loadset(filename);
        % Interpolation of predefined channels
        EEG = interpChan(EEG,analysis,iSubj,trig.triggers{iTrig,1});
        % Select all events marked as iTrig
        EEG2 = pop_selectevent(EEG,'type',trig.triggers{iTrig,1},'latency','-10 <= 10','deleteevents','off','deleteepochs','on');
        clear EEG;
        % display number of events for the current trigger and the current subject before rejection
        disp(' ');
        disp([':: Subject ' num2str(iSubj, '%0.2d') ...
              ' Trigger ' num2str(trig.triggers{iTrig,1}) ' ' num2str(size(EEG2.data,3)) ' Events original']);
        
        switch 1
          case analysis.rejmode == 0
            % NO REJECTION
            disp(':: ARTIFACT REJECTION IS SWITCHED OFF!');
          case analysis.rejmode == 1 | analysis.rejmode == 2
            % DELTA REJECTION (with or without eye) go through all channels in all trials
            % and count events that will be rejected according to delta
            % criterion
            for iChannel = 1:size(EEG2.data,1)
                thresh = 0;
                nRej = 0;
                for iTrial = 1:size(EEG2.data,3)
                    thresh = abs(max(EEG2.data(iChannel,:,iTrial))- ...
                                 min(EEG2.data(iChannel,:,iTrial)));

                    
                    if thresh > analysis.sortthresh
                        nRej = nRej + 1;
                    end
                end
                % If more than the given proportion threshold threshold of
                % rejected epoches on this channel is exceeded -> report
                % channel in rejFile
                if (nRej/size(EEG2.data,3)) > analysis.chanMaxRej
                    fprintf(rejFile,'%s\n',[':: Subject ' num2str(iSubj,'% 0.2d') ...
                                        ' Trigger ' num2str(trig.triggers{iTrig,1}) ...
                                        ' channel ' num2str(iChannel) ' (' EEG2.chanlocs(iChannel).labels ') '...
                                        'rejected Epochs: ' num2str((nRej/size(EEG2.data,3))) ...
                                        ' %']);
                end
            end
            % reject epochs with minimal absolut amplitude differences (due to technical
            % problems)
            if analysis.rejFlatepochs
                epochOut = zeros(1,size(EEG2.data,3));
                minChange = 100;
                % find events with minimal amplitude change
                for iChannel = 1:size(EEG2.data,1)
                    for iTrial = 1:size(EEG2.data,3)
                        amplChange = abs(max(EEG2.data(iChannel,:,iTrial)) - min(EEG2.data(iChannel,:,iTrial)));
                        if amplChange < analysis.flatthresh
                            epochOut(iTrial) = 1;                                
                        end
                        minChange = min(minChange,amplChange);
                    end
                end
                % count number of events with minimal amplitude change
                nFlat = numel(find(epochOut==1));
                if nFlat
                    disp([':: Subject ' num2str(iSubj, '%0.2d') ...
                          ' Trigger ' num2str(trig.triggers{iTrig,1}) ' ' num2str(nFlat) ' Events flat']);
                    disp(':: These will be rejected');
                    % reject events with minimal amplitude change                       
                    EEG2.reject.rejmanual = epochOut;
                    % report to file
                    fprintf(rejFile,'%s\n',['           theoretically ' num2str(eventCount(iTrig)) ' trials, practically ' num2str(size(EEG2.data,3)) ' trials']);
                    if analysis.rejFlatepochs && nFlat
                        fprintf(rejFile,'%s\n',['           ' num2str(nFlat) ' trials rejected for being flat']);
                    end
                    disp([':: Subject ' num2str(iSubj, '%0.2d') ' Trigger ' num2str(trig.triggers{iTrig,1}) ' wrong number of trials']);
                    
                else
                    disp([':: Subject ' num2str(iSubj, '%0.2d') ...
                          ' Trigger ' num2str(trig.triggers{iTrig,1}) [' minimum ' ...
                                        'amplitude change: '] num2str(minChange) ' microV']);
                end
            end 
            
            % mark trials for rejection
            EEG2 = eeg_rejdelta(EEG2,'thresh',analysis.sortthresh,'chans',setdiff(1:size(EEG2.data,1),exclVector));
          case analysis.rejmode == 3
            % REJECTION BASED ON PREDEFINED TRIAL INDICES 
            
            % if loop is entered the first time load the rejection index matrix
            if iTrig == 1
                rejEpoch = rejIndex(paths);
            end
          otherwise
            error([':: Invalid Option for rejection mode: ' num2str(analysis.rejmode)]);
        end
        
        % Superpose all marked rejections to reject.rejglobal
        EEG2 = eeg_rejsuperpose(EEG2, 1, 1, 1, 1, 1, 1, 1, 1);
        % get vector of marked trials and store them for each participant and each trigger
        switch 1
          case analysis.rejmode == 3
            % assign rejection index vector of current subject and current trigger to
            % EEG2.reject.rejglobal vector
            EEG2.reject.rejglobal = rejEpoch{iSubj,iTrig};
          case analysis.rejmode == 1 | analysis.rejmode == 2
            % get vector of marked trials and store them for each participant and each trigger
            rejEpoch{1,iTrig} = EEG2.reject.rejglobal;
        end            
        
        %save good and bad trials
        EEG3 = pop_rejepoch(EEG2,EEG2.reject.rejglobal,0);
        EEG3 = pop_saveset(EEG3,[...
            paths.resFileSubSpec num2str(iSubj, '%0.2d')... 
            paths.resFileTrigSpec num2str(trig.triggers{iTrig,1},'%0.2d') 'good' '.set'],subjDir);

        %% get ERPs for current Subject %%%%%%%%%%%%%%%%%%%%%
        subErp(1,iTrig,:,:) = (squeeze(mean(EEG3.data,3)))';
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            

        disp([':: Subject ' num2str(iSubj, '%0.2d') ' Trigger ' num2str(trig.triggers{iTrig,1}) ' after correction: ' num2str(size(EEG3.data,3)) ' Events averaged']);
        % save number of usable trials
        corrTrials(1,iTrig) = size(EEG3.data,3);
        
        % only 1 usable epoch, probably created artificially
        if corrTrials(1,iTrig) == 1 
            corrTrials(1,iTrig) = 0;
        end
        
        disp(' ');
        
        disp(':: Extracting and saving bad trials.');
        EEG3 = pop_rejepoch(EEG2,mod(EEG2.reject.rejglobal+1,2),0);
        EEG3 = pop_saveset(EEG3,[...
            paths.resFileSubSpec num2str(iSubj, '%0.2d')... 
            paths.resFileTrigSpec num2str(trig.triggers{iTrig,1},'%0.2d') 'bad' '.set'],subjDir);
        % save chanlocs
        if ~exist([paths.resDirAll 'chanlocs.mat'],'file')            
            disp(' ');
            disp(':: Saving channel location information temporarily'); 
            chanlocs = EEG2.chanlocs;
            save([paths.resDirAll 'chanlocs.mat'],'chanlocs');
            clear chanlocs;
        end
        
        % clear EEG structures
        clear EEG2;
        clear EEG3;
        
    end % if current file exists
    close all;
end % trigger loop
% close file
fclose(rejFile);

%% Get Event Indices of associated Triggers
maxIndex = max(cell2mat(trig.triggers(:,6)));
trialNum = cell(maxIndex,1);
for iInd = 1:maxIndex
    for iTrig = 1:size(trig.triggers,1)
        if trig.triggers{iTrig,6} == iInd & corrTrials(1,iTrig)
           trialNum{iInd,1} = [trialNum{iInd,1} corrTrials(1,iTrig)];
        end
    end
end

% get list of  file names in current dir
subjFiles = dir([subjDir '*.set']);

% initialize cell array of randomly chosen trial indices
subTrialInd = cell(size(trig.triggers,1),1);
for iTrig = 1:size(trig.triggers,1)
    % get current file
    filename = [subjDir paths.resFileSubSpec num2str(iSubj, '%0.2d')...
                paths.resFileTrigSpec num2str(trig.triggers{iTrig,1},'%0.2d') 'good.set'];
    if exist(filename,'file')
        % load current file
        EEG = pop_loadset(filename);
        % get number of trials to average
        if ~isempty(trig.triggers{iTrig,6})
            % if the current trial is associated with other trials, get the
            % number of trials for the current trial (leading to equal
            % [minimal] number of trials for all the trials in the
            % associated trigger pool)

            % randomly choose currTrialNum times epochs of the current
            % trial
            availTrials = 1:size(EEG.data,3);
            trialInd = [ ];
            
            for iTrial = 1:min(trialNum{trig.triggers{iTrig,6}})
                % randomly pick one epoch
                if ~isempty(availTrials)
                    currTrial = randi(length(availTrials));
                    trialInd = [trialInd  availTrials(currTrial)];
                    % remove the trial from available trials
                    availTrials(currTrial) = [ ];       
                end
            end
            % pick the chosen epochs from EEG structure
            if ~isempty(trialInd)
                EEG.data = EEG.data(:,:,trialInd);
            else
                EEG.data(:,:,:) = 0;
            end
                
        else 
            % if the current trigger is not associated with other triggers,
            % take all epochs
            trialInd = [ ];
        end
        % save current trial indices
        subTrialInd{iTrig,1} = trialInd;
        
        %% get ERPs for current Subject %%%%%%%%%%%%%%%%%%%%%
        subErpEqual(1,iTrig,:,:) = (squeeze(mean(EEG.data,3)))';
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        clear EEG;
    end
end

% delete current subject folder
if analysis.clearFolders
    rmdir(subjDir,'s')    
end
