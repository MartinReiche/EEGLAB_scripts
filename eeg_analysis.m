% Metascript to adjust analysis parameters and run through all analysis
% steps using EEGLAB
%
% possible rejection candidates: ??
%
% Needs: 
% eeg_emcp - Andreas Widmann (http://github.com/widmann)
% firfilt 1.5 - Andreas Widmann (http://github.com/widmann)
% eeg_rejdelta - Andreas Widmann (http://github.com/widmann)
% RMAOV1 - Trujillo-Ortiz, A., R. Hernandez-Walls and
%      R.A. Trujillo-Perez. (2004). RMAOV1:One-way repeated measures ANOVA. A
%      MATLAB file. [WWW document]. URL
%      http://www.mathworks.com/matlabcentral/
%      fileexchange/loadFile.do?objectId=5576
%
% Files you need to configure:
% - config (central configuration file)
% - triggerlabels.m (configuratiuon of trigger labels)
% - func/retrigConf.m (configuration for retriggering)
% - func/changeTrig.m (define trigger range for omissions [for exclusions
%   around those events])
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

% Please define all parameters in config.m
function eeg_analysis(taskType,subjects,analysisIn,filtParIn,trigIn,batchMode)

    tStart = tic; % start measuring time for Analysis
    
    %% Get Parameters
    if nargin > 2
        % Get analysis parameters from input arguments [batch mode]
        analysis = analysisIn;
        % save subjects vector in analysis struct from input arguments [batch mode]
        analysis.subjects = subjects;
        % Get filter parameters from input arguments [batch mode]
        filtPar = filtParIn;
        % Load trigger Matrix for current task from input arguments [batch mode]
        trig = trigIn;
        analysis.batchMode = 1;
        analysis.preprocess = 1;
    else
        % Get analysis parameters
        analysis = config('parameters','task',taskType);
        % save subjects vector in analysis struct
        analysis.subjects = subjects;
        % Get filter parameters
        filtPar = config('filter');
        % Load trigger Matrix for current task
        trig = config('triggers','task', taskType);
        analysis.jobIndex = 0;
        analysis.batchMode = 0;
        analysisIn = 0;
        filtParIn = 0;
        trigIn = 0;
        batchMode = 0;
    end    

    %% Add Paths
    % get Paths
    paths = config('Path','task',taskType,'analysis',analysis,'filt',filtPar);
    % add all relevant paths to the current workspace
    addpath(paths.funcDir, paths.libDir, paths.eeglabDir, paths.stimFuncDir);
    addpath([paths.eeglabDir 'functions/popfunc/']);
    addpath([paths.eeglabDir 'functions/guifunc/']);
    addpath([paths.eeglabDir 'functions/adminfunc/']);
    addpath([paths.eeglabDir 'functions/sigprocfunc/']);
    addpath([paths.eeglabDir 'external/biosig-partial/t200_FileAccess/']);
    addpath([paths.eeglabDir 'external/biosig-partial/t250_ArtifactPreProcessingQualityControl/']);
    addpath([paths.eeglabDir 'plugins/bva-io1.58/']);
    addpath([paths.libDir 'firfilt-1.5.3/']);
    addpath([paths.libDir 'sphspline0.2/']);

    % check availabel raw data and configure raw data paths
    paths = checkRawData(paths,subjects,taskType);
    
    %% run preprocessing routine
    if analysis.preprocess
        % initialize cell array for subject erps
        subErp = cell(numel(subjects),1);
        subErpEqual = cell(numel(subjects),1);
        subTrialInd = cell(numel(subjects),1);
        corrTrials = cell(numel(subjects),1);
        numEvent = cell(numel(subjects),1);
        trigNum = cell(numel(subjects),1);
        % open independant matlab labs (local configuration [4 cores])
       
        if analysis.parallel
            % open matlabpool according to specified parallel computing profile)
            matlabpool(analysis.core)
        end
        
        % over all subjects
        parfor iSubj = 1:size(subjects,2)

            %% Get Parameters
            if batchMode
                % Get analysis parameters from input arguments [batch mode]
                analysis = analysisIn;
                % save subjects vector in analysis struct from input arguments [batch mode]
                analysis.subjects = subjects;
                % Get filter parameters from input arguments [batch mode]
                filtPar = filtParIn;
                % Load trigger Matrix for current task from input arguments [batch mode]
                trig = trigIn;
                analysis.batchMode = 1;
            else
                % Get analysis parameters
                analysis = config('parameters','task',taskType);
                % save subjects vector in analysis struct
                analysis.subjects = subjects;
                % Get filter parameters
                filtPar = config('filter');
                % Load trigger Matrix for current task
                trig = config('triggers','task', taskType);
                analysis.jobIndex = 0;
                analysis.batchMode = 0;
            end    

            
            %% Add Paths
            % get Paths
            paths = config('Path','task',taskType,'analysis',analysis,'filt',filtPar);

            addpath(paths.funcDir, paths.libDir, paths.eeglabDir, paths.stimFuncDir);
            addpath([paths.eeglabDir 'functions/popfunc/']);
            addpath([paths.eeglabDir 'functions/guifunc/']);
            addpath([paths.eeglabDir 'functions/adminfunc/']);
            addpath([paths.eeglabDir 'functions/sigprocfunc/']);
            addpath([paths.eeglabDir 'external/biosig-partial/t200_FileAccess/']);
            addpath([paths.eeglabDir 'external/biosig-partial/t250_ArtifactPreProcessingQualityControl/']);
            addpath([paths.libDir 'firfilt-1.5.3/']);
            addpath([paths.libDir 'sphspline0.2/']);

            % check availabel raw data and configure raw data paths
            paths = checkRawData(paths,subjects,taskType);

            % % start measuring time for Analysis
            % tStart = tic;

            % save subjects vector in analysis struct
            analysis.subjects = subjects;
            % prepare matrix for counting of available events for each subjects
            % (lines represent blocks, colums 1-20 are triggers 101-512)
            eventCount = zeros(size(trig.triggers,1),analysis.nBlocks);
            % prepare result folder for subject
            paths = prepSubDir(paths,subjects,iSubj);

            % get condition order
            condOrder = preporder(subjects(iSubj),taskType,0);
            counter = 0;
            for iFile = 1:analysis.nBlocks
                counter = counter + 1;
                % load raw data for current file and stimulation parameters
                [EEG, block,analysis, paths] = loadRawData(paths, subjects(iSubj),taskType,iFile,analysis,counter);
                % check parameters (sampling rate, triggers etc)
                switch analysis.rawFormat
                  case 'biosemi'
                    [EEG, eventExcp] = checkFileBiosig(EEG,subjects(iSubj),iFile,block,taskType,analysis,trig);
                  case 'brainvision'
                    [EEG, eventExcp] = checkFileBv(EEG,subjects(iSubj),iFile,block,taskType,analysis,trig);
                end
                % Retriggering and systematically exclude events from analysis
                EEG = change_trig(EEG,analysis,trig,iFile,condOrder,taskType); 
                % bipolarize eye channels
                EEG = bipolarize(EEG,analysis);
                % adjust electrode positions
                EEG = chanLoc(EEG,paths);
                % Run Channel interpolation
                EEG = interpChan(EEG,analysis,subjects(iSubj),iFile);
                % filter the block depending on rejection method:
                switch analysis.rejmode
                    % No eye correction
                  case {0,1,3,4}
                    % filter Data
                    filtPar.eye = 0; % no eye correction
                    EEG = fir_filter(EEG,analysis,filtPar);
                    
                    % safe block data after segemtation and baseline correction
                    eventCount = segmentation(EEG,trig,analysis,paths,eventCount, ...
                                              iFile,subjects(iSubj),condOrder);
                    % eye correction
                  case 2
                    % filter Data
                    filtPar.eye = 1; % before eye correction
                    EEG = fir_filter(EEG,analysis,filtPar);
                    % save EEG block
                    pop_saveset(EEG,[paths.resFileSubSpec num2str(subjects(iSubj),'%0.2d') ...
                                     paths.resFileBlockSpec num2str(iFile, '%0.2d') '.set'], ...
                                paths.resDir);
                  otherwise
                    error([':: Invalid option for rejection method: ' num2str(analysis.rejmode)]);
                end
            end
                
            switch analysis.rejmode
                % no eye correction
              case {0,1,3}
                % merge files of same trigger for current subject and save
                merge_data(subjects(iSubj),paths,trig,eventCount,condOrder);
                % eye correction
              case 2
                % merge all files for current subject, perform eye correction, apply post
                % filter and save
                merge_all(subjects(iSubj),paths,filtPar,analysis);
                % safe block data after trigger resegemtation and baseline correction
                segmentation(EEG,trig,analysis,paths,eventCount, ...
                             iFile,subjects(iSubj),condOrder);
            end
            % artifact rejection and averaging
            [subErp{iSubj,1},...
             subErpEqual{iSubj,1},...
             subTrialInd{iSubj,1},...
             corrTrials{iSubj,1},...
             numEvent{iSubj,1},...
             rejEpoch(iSubj,:),...
             trialNum{iSubj,1}] = eeg_rejection(subjects(iSubj),paths,trig,analysis);
        end % parfor all subjects
        clear EEG

        if analysis.parallel
            % close the pool
            matlabpool close; 
        end
        
        % combine subject data, calculate amount of usable data and save
        % erpAll
        tEnd = toc(tStart);
        save_erp(subErp,subErpEqual,subTrialInd,corrTrials,numEvent,rejEpoch,trialNum,subjects,paths,trig,tStart,tEnd,analysis,filtPar);

        disp(' ');
        disp([':: Duration of Analysis: ' num2str(tEnd/60) ' minutes.']);
        % clear current trigger configuration (will be loaded from erp file for
        % further steps)
        clear trig;
    end % preprocess

    if ~batchMode
        %% Calculation of difference
        % calculate difference waves
        [erpAll, restoredConf, chanlocs, trig] = calc_diff(subjects,paths,analysis,taskType);
        % Get Plot Parameters
        [plotPar] = config('Plot','task',taskType);
        % Select Channels to plot
        plot_erp(erpAll,chanlocs,plotPar,trig,analysis,paths,taskType,restoredConf);
    end

    
    