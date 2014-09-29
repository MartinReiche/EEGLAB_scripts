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

function [subErp,subErpEqual,subTrialInd,corrTrials,eventCount,rejEpoch,trialNum,rejLog] = eeg_rejection(iSubj,paths,trig,analysis)
% get current dir
subjDir = [paths.resDirAll paths.resSubFolderPrefix num2str(iSubj, '%0.2d') '/'];
% get list of  file names in current dir
subjFiles = dir([subjDir '*set']);
% load one (the 1st) of the .set files for the current subject
EEG = pop_loadset(subjFiles(1).name,subjDir);
% initialize erp matrix for curren subject
subErp = zeros(1,size(trig.triggers,1),size(EEG.data,2),size(EEG.data,1));
subErpEqual = zeros(1,size(trig.triggers,1),size(EEG.data,2),size(EEG.data,1));
% initialize cell arrays for vector of rejected epoch indices per trial
rejEpoch = cell(1,size(trig.triggers,1));
% get chanlocs
chanlocs = EEG.chanlocs;
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

% initialize matrix for number of usable data per trigger
corrTrials = zeros(1,size(trig.triggers,1));

%% Start Rejection for all files
% initialize cell for logging rejection specific information
rejLog.chanReport = {};
rejLog.flat = {};
rejLog.minChange = zeros(size(trig.triggers,1),1);
rejLog.maxChange = rejLog.minChange;

for iTrig = 1:size(trig.triggers,1)
    % convert cell of trigger names to concatenated trigger string in
    % case several triggers are selected at once
    if iscell(trig.triggers{iTrig,1})
        currTrig = num2str(cell2mat(trig.triggers{iTrig,1}),'%0.2d');
    else
        currTrig = num2str(trig.triggers{iTrig,1},'%0.2d');
    end
    
    % get current file
    filename = [subjDir paths.resFileSubSpec num2str(iSubj, '%0.2d')...
                paths.resFileTrigSpec currTrig 'all.set'];
    % if the current epoch file exists
    if exist(filename,'file')
        EEG = pop_loadset(filename);
        % Interpolation of predefined channels
        EEG = interpChan(EEG,analysis,iSubj,currTrig);
        % Select all events marked as iTrig
        EEG2 = pop_selectevent(EEG,'type',trig.triggers{iTrig,1},'latency','-10 <= 10','deleteevents','off','deleteepochs','on');
        
        clear EEG;
        % display number of events for the current trigger and the current subject before rejection
        disp(' ');
        disp([':: Subject ' num2str(iSubj, '%0.2d') ...
              ' Trigger ' currTrig ' ' num2str(size(EEG2.data,3)) ' Events original']);
        
        switch 1
          case analysis.rejmode == 0
            % NO REJECTION
            disp(':: ARTIFACT REJECTION IS SWITCHED OFF!');
          case ismember(analysis.rejmode,[1 2])
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
                    % save for later report in rejection log file [handled
                    % by save_erp.m]
                    rejLog.chanReport = [rejLog.chanReport;{iSubj currTrig iChannel ...
                                        EEG2.chanlocs(iChannel).labels (nRej/size(EEG2.data,3))}];
                end
            end

            % get minimal and maximal amplitude values aover all channels
            % and trials for current trigger and determine flat epoches of
            % current trigger 
            epochOut = zeros(1,size(EEG2.data,3));
            amplChange = zeros(size(EEG2.data,1),size(EEG2.data,3));
            % find events with minimal amplitude change
            for iChannel = 1:size(EEG2.data,1)
                for iTrial = 1:size(EEG2.data,3)
                    amplChange(iChannel,iTrial) = abs(max(EEG2.data(iChannel,:,iTrial)) - min(EEG2.data(iChannel,:,iTrial)));
                    if amplChange(iChannel,iTrial) < analysis.flatthresh
                        epochOut(iTrial) = 1;                                
                    end
                end
            end
            rejLog.minChange(iTrig) = min(min(amplChange));
            rejLog.maxChange(iTrig) = max(max(amplChange));
            % count number of events with minimal amplitude change
            nFlat = numel(find(epochOut == 1));
                        
            % reject epochs with minimal absolut amplitude differences (due to technical
            % problems)
            if analysis.rejFlatepochs && nFlat
                disp([':: Subject ' num2str(iSubj, '%0.2d') ...
                      ' Trigger ' currTrig ' ' num2str(nFlat) ' Events flat']);
                disp(':: These will be rejected');
                % reject events with minimal amplitude change                       
                EEG2.reject.rejmanual = epochOut;
                % store number of trial which were rejected for being
                % flat for later report in rejection log [handled by save_erp.m]
                rejLog.flat = [rejLog.flat;{iSubj iTrig nFlat eventCount(iTrig) size(EEG2.data,3)}];
                disp([':: Subject ' num2str(iSubj, '%0.2d') ' Trigger ' currTrig ' wrong number of trials']);
            else
                disp([':: Subject ' num2str(iSubj, '%0.2d') ...
                      ' Trigger ' currTrig [' minimum ' ...
                                    'amplitude change: '] num2str(rejLog.minChange(iTrig)) ' micro Volts']);
                disp(' ');
            end

            
            % mark trials for rejection
            EEG2 = eeg_rejdelta(EEG2,'thresh',analysis.sortthresh,'chans',1:size(EEG2.data,1));
          case analysis.rejmode == 3
            % REJECTION BASED ON PREDEFINED TRIAL INDICES 
            
            % if loop is entered the first time load the rejection index matrix
            if iTrig == 1
                rejEpoch = rejIndex(paths);
            end
          case ismember(analysis.rejmode,[4 5])
            % SORTED AVERAGING WITH AND WITHOUT EYE CORRECTION
                        
            % BASELINE CORRECTION            
            disp(':: Performing baseline correction'); 
            % get ms range of baseline window relative to beginning of epoch
            baseMS(1) = analysis.baseWin(1) - analysis.erpWin(1);
            baseMS(2) = analysis.baseWin(2) - analysis.erpWin(1);
            % correct for 0
            baseMS(baseMS == 0) = 1;
            % get the point range of the baseline window
            pointrange = ceil(baseMS(1)*analysis.sampRate/1000):floor(baseMS(2)*analysis.sampRate/1000);
            % separately for all subjects and all channels
            for iChan = 1:size(EEG2.data,1)
                % get the mean value in the baseline window of the current subject on the
                % current channel and subtract it from the whole data range of the
                % current subject on the current channel
                tmpmean = mean(double(EEG2.data(iChan,pointrange,:)),2);
                % erpAll(iSubj,:,:,iChan) = erpAll(iSubj,:,:,iChan) - repmat(tmpmean, [1 1 size(erpAll,3) 1]);
                EEG2.data(iChan,:,:) = EEG2.data(iChan,:,:) - repmat(tmpmean, [1 size(EEG2.data,2) 1]);
                % mean(double(EEG2.data(iChan,pointrange,:)),2)
            end
            
            disp(':: Performing sorted averaging'); 
            % SORTED AVERAGING (with or without eye movement correction) 
            % mark trials for rejection
            EEG2 = eeg_rejsortavg(EEG2);

            % get minimal and maximal amplitude values aover all channels
            % and trials for current trigger and determine flat epoches of
            % current trigger 
            epochOut = zeros(1,size(EEG2.data,3));
            amplChange = zeros(size(EEG2.data,1),size(EEG2.data,3));
            % find events with minimal amplitude change
            for iChannel = 1:size(EEG2.data,1)
                for iTrial = 1:size(EEG2.data,3)
                    amplChange(iChannel,iTrial) = abs(max(EEG2.data(iChannel,:,iTrial)) - min(EEG2.data(iChannel,:,iTrial)));
                    if amplChange(iChannel,iTrial) < analysis.flatthresh
                        epochOut(iTrial) = 1;                                
                    end
                end
            end
            rejLog.minChange(iTrig) = min(min(amplChange));
            rejLog.maxChange(iTrig) = max(min(amplChange));
            % count number of events with minimal amplitude change
            nFlat = numel(find(epochOut == 1));
                        
            % reject epochs with minimal absolut amplitude differences (due to technical
            % problems)
            if analysis.rejFlatepochs && nFlat
                disp([':: Subject ' num2str(iSubj, '%0.2d') ...
                      ' Trigger ' currTrig ' ' num2str(nFlat) ' Events flat']);
                disp(':: These will be rejected');
                % reject events with minimal amplitude change                       
                EEG2.reject.rejmanual = epochOut == 1 | EEG.reject.rejmanual == 1;
                % store number of trial which were rejected for being
                % flat for later report in rejection log [handled by save_erp.m]
                rejLog.flat = [rejLog.flat;{iSubj iTrig nFlat eventCount(iTrig) size(EEG2.data,3)}];
                disp([':: Subject ' num2str(iSubj, '%0.2d') ' Trigger ' currTrig ' wrong number of trials']);
            else
                disp([':: Subject ' num2str(iSubj, '%0.2d') ...
                      ' Trigger ' currTrig [' minimum ' ...
                                    'amplitude change: '] num2str(rejLog.minChange(iTrig)) ' microV']);
                disp(' ');
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
          case ismember(analysis.rejmode,[1 2 4 5])
            % get vector of marked trials and store them for each participant and each trigger
            rejEpoch{1,iTrig} = EEG2.reject.rejglobal;
        end            
        
        %save good and bad trials
        EEG3 = pop_rejepoch(EEG2,EEG2.reject.rejglobal,0);
        EEG3 = pop_saveset(EEG3,[...
            paths.resFileSubSpec num2str(iSubj, '%0.2d')... 
            paths.resFileTrigSpec currTrig 'good' '.set'],subjDir);

        %% get ERPs for current Subject %%%%%%%%%%%%%%%%%%%%%
        subErp(1,iTrig,:,:) = (squeeze(mean(EEG3.data,3)))';
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
        disp(' ');
        disp([':: Subject ' num2str(iSubj, '%0.2d') ' Trigger ' currTrig ' after correction: ' num2str(size(EEG3.data,3)) ' Events averaged']);
        % save number of usable trials
        corrTrials(1,iTrig) = size(EEG3.data,3);
        
        % only 1 usable epoch, probably created artificially
        if corrTrials(1,iTrig) == 1 
            corrTrials(1,iTrig) = 0;
        end
        disp(':: Extracting and saving bad trials.');
        disp(' ');
        EEG3 = pop_rejepoch(EEG2,mod(EEG2.reject.rejglobal+1,2),0);
        EEG3 = pop_saveset(EEG3,[...
            paths.resFileSubSpec num2str(iSubj, '%0.2d')... 
            paths.resFileTrigSpec currTrig 'bad' '.set'],subjDir);
        % save chanlocs
        if ~exist([paths.chanlocs 'chanlocs.mat'],'file')            
            disp(' ');
            disp(':: Saving channel location information temporarily'); 
            chanlocs = EEG2.chanlocs;
            save([paths.chanlocs 'chanlocs.mat'],'chanlocs');
            clear chanlocs;
        end
        
        % clear EEG structures
        clear EEG2;
        clear EEG3;
        
    end % if current file exists
    close all;
end % trigger loop

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

    % convert cell of trigger names to concatenated trigger string in
    % case several triggers are selected at once
    if iscell(trig.triggers{iTrig,1})
        currTrig = num2str(cell2mat(trig.triggers{iTrig,1}),'%0.2d');
    else
        currTrig = num2str(trig.triggers{iTrig,1},'%0.2d');
    end
    % get current file
    filename = [subjDir paths.resFileSubSpec num2str(iSubj, '%0.2d')...
                paths.resFileTrigSpec currTrig 'good.set'];
    
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
