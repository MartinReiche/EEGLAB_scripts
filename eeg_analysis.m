% Metascript to adjust analysis parameters and run through all analysis
% steps using EEGLAB
%
% USAGE: 
% in serial mode:       eeg_analysis(taskType,[numbers of subjects])
% in parallel mode:     eeg_analysis(taskType,[numbers of subjects],'parallelMethod')
% 
%  'parallelMethod':    'submitJob' - submit prprocessing job to cluster
%                       'getJob' - retrieve preprocessed job data from cluster
%                       'plot' - omit preprocessing and plot processed data
%
% NEEDS:
% 
% EEGLAB 13    - A Delorme & S Makeig (2004) EEGLAB: an open source toolbox 
%                for analysis of single-trial EEG dynamics. 
%                Journal of Neuroscience Methods 134:9-21 
% bv-io        - can be installed via EEGLAB extension manager 
% eeg_emcp     - Andreas Widmann (http://github.com/widmann)
% firfilt 1.6  - Andreas Widmann (http://github.com/widmann)
% eeg_rejdelta - Andreas Widmann (http://github.com/widmann)
% RMAOV1       - Trujillo-Ortiz, A., R. Hernandez-Walls and
%                R.A. Trujillo-Perez. (2004). RMAOV1:One-way repeated measures ANOVA. A
%                MATLAB file. URL: http://www.mathworks.com/matlabcentral/
%                fileexchange/loadFile.do?objectId=5576
%
% Files you need to configure:
% - config (central configuration file)
% - triggerlabels.m (configuratiuon of trigger labels)
% - channelInterp.m (configuration of to-be-interpolated channels for given subjects)
%
% Copyright (c) 2014 Martin Reiche, Carl-von-Ossietzky-University Oldenburg
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
function [sched,job] = eeg_analysis(taskType,subjects,parallelMode)

    % Get analysis parameters
    analysis = config('parameters','task',taskType);
    % save subjects vector in analysis struct
    analysis.subjects = subjects;
    % Get filter parameters
    filtPar = config('filter');
    % Load trigger Matrix for current task
    trig = config('triggers','task', taskType);
    % check input
    if analysis.parallel && nargin < 3
        error(':: Parallel Mode is required');
    elseif nargin > 2
        analysis.parallel = 1;
    end
    %% Add Paths
    % get Paths
    paths = config('Path','task',taskType,'analysis',analysis,'filt',filtPar);
    % add local paths
    addpath(paths.funcDir, paths.local.libDir, paths.local.eeglabDir, paths.local.stimFuncDir);
    addpath([paths.local.eeglabDir 'functions/popfunc/']);
    addpath([paths.local.eeglabDir 'functions/guifunc/']);
    addpath([paths.local.eeglabDir 'functions/adminfunc/']);
    addpath([paths.local.eeglabDir 'functions/sigprocfunc/']);
    addpath([paths.local.eeglabDir 'plugins/bva-io1.5.12/']);
    addpath([paths.local.eeglabDir 'plugins/Biosig2.88/']);
    addpath([paths.local.eeglabDir 'plugins/Biosig2.88/biosig/doc']);
    addpath([paths.local.eeglabDir 'plugins/Biosig2.88/biosig/t250_ArtifactPreProcessingQualityControl']);
    addpath([paths.local.eeglabDir 'plugins/Biosig2.88/biosig/t200_FileAccess']);
    addpath([paths.local.eeglabDir 'plugins/firfilt1.6/']);
    addpath([paths.local.libDir 'sphspline0.2/']);
    
    %% run preprocessing routine
    if analysis.preprocess
        % evaluate runmode (parallel or serial)
        if ~analysis.parallel && nargin < 3
            %% serial processing on local machine
            % check availabel raw data and configure raw data paths
            paths = checkRawData(paths,subjects,taskType);
            % initialize cell array for subject erps
            subErp = cell(numel(subjects),1);
            subErpEqual = cell(numel(subjects),1);
            subTrialInd = cell(numel(subjects),1);
            corrTrials = cell(numel(subjects),1);
            numEvent = cell(numel(subjects),1);
            trigNum = cell(numel(subjects),1);
            % over all subjects
            dur.Start = datestr(now,'ddd mmm DD HH:MM:SS YYYY');
            for iSubj = 1:size(subjects,2)
                dur.task(iSubj).Start = datestr(now,'ddd mmm DD HH:MM:SS YYYY');
                subData = preprocess(taskType,subjects(iSubj),analysis,filtPar,trig,paths);
                % combine subject specific results of preprocessing
                subErp{iSubj,1} = subData.subErp;
                subErpEqual{iSubj,1} = subData.subErpEqual;
                subTrialInd{iSubj,1} = subData.subTrialInd;
                corrTrials{iSubj,1} = subData.corrTrials;
                numEvent{iSubj,1} = subData.numEvent;
                rejEpoch(iSubj,:) = subData.rejEpoch;
                trialNum{iSubj,1} = subData.trialNum;                 
                rejLog(iSubj) = subData.rejLog;
                dur.task(iSubj).End = datestr(now,'ddd mmm DD HH:MM:SS YYYY');
            end 
            dur.End = datestr(now,'ddd mmm DD HH:MM:SS YYYY');
            % retrieve configuration settings
            analysis = subData.analysis;
            paths = subData.paths;
            trig = subData.trig;
            filtPar = subData.filtPar;
            % save preprocessed data
            save_erp(subErp,subErpEqual,subTrialInd,corrTrials,numEvent,rejEpoch,trialNum,subjects,paths,trig,dur,analysis,filtPar,rejLog);
            job = [];
            sched = [];
          else
            %% parallel processing
            
            % Create scheduler %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % all jobs defined below, will be scheduled by this cluster configuration
            sched = findResource('scheduler','configuration',analysis.core);
            % set resources for scheduler
            if strcmp(paths.cluster,'remote')
                set(sched, 'SubmitFcn', cat(1,sched.SubmitFcn,'runtime','0:30:0','memory','6G'));
            end
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            switch parallelMode
              case 'submitJob'
                % create cell array of input argument for each task
                inArgs = cell(1,size(subjects,2));
                for iTask = 1:size(subjects,2)
                    inArgs{iTask} = {taskType subjects(iTask),analysis,filtPar,trig,paths};
                end
                % get list of file dependencies (in main folder)
                files = dir([pwd '/*.m']);
                fd1 = cell(1,size(files,1));
                for iFile = 1:size(files,1)
                    fd1{iFile} = [pwd '/' files(iFile).name];
                end
                % and in functions folder
                files = dir([pwd '/func/*.m']);
                fd2 = cell(1,size(files,1));
                for iFile = 1:size(files,1)
                    fd2{iFile} = [pwd '/func/' files(iFile).name];
                end
                % combine file dependencies from func and main folder
                fd = [fd1 fd2];
                % create job object 
                job = createJob(sched);
                % set file dependencies 
                set(job,'FileDependencies',fd);
                task = createTask(job, @preprocess, 1 ,inArgs);

                % add paths dependencies
                pd = {paths.rawDir; paths.resDir; paths.behavDir; ...
                      paths.eeglabDir; paths.libDir; paths.elecSetup;...
                      paths.stimFuncDir;
                      [paths.eeglabDir 'functions/popfunc/'];...
                      [paths.eeglabDir 'functions/guifunc/'];...
                      [paths.eeglabDir 'functions/adminfunc/'];...
                      [paths.eeglabDir 'functions/sigprocfunc/'];...
                      [paths.eeglabDir 'plugins/bva-io1.5.12/'];...
                      [paths.eeglabDir 'plugins/plugins/Biosig2.88/'];...
                      [paths.eeglabDir 'plugins/Biosig2.88/biosig/doc'];...
                      [paths.eeglabDir 'plugins/Biosig2.88/biosig/t250_ArtifactPreProcessingQualityControl'];...
                      [paths.eeglabDir 'plugins/Biosig2.88/biosig/t200_FileAccess'];...
                      [paths.eeglabDir 'plugins/firfilt1.6/'];...
                      [paths.libDir 'sphspline0.2/']};                
                
                set(job,'PathDependencies',pd);

                % submit the job to the scheduler
                submit(job);
                
              case 'getJob'
                job = get(sched,'Jobs');
                % specify Job ID
                if size(job,1) > 1
                    job
                    disp(':: Specify job ID')
                    id = 0;
                    while ~ismember(id,1:size(job,1))
                        id = input('>> ');
                        if ~ismember(id,1:size(job,1))
                           disp([':: ' num2str(id) ' is not a valid job ID']);  
                        end
                    end
                    job = job(id);
                elseif size(job,1) < 1
                    error(':: There are no jobs for the current scheduler');
                end
                % get Tasks
                task = get(job,'Tasks');
                % check task state
                allFin = [ ];
                for iTask = 1:size(task,1)
                    if strcmp(task(iTask).State,'finished')
                        allFin = [allFin 1];
                    else
                        allFin = [allFin 0];
                    end
                end
                % gather data
                if all(allFin) && ~isempty(allFin)
                    disp(':: All tasks are finished');
                    task
                    
                    % initialize cell array for subject erps
                    subErp = cell(numel(subjects),1);
                    subErpEqual = cell(numel(subjects),1);
                    subTrialInd = cell(numel(subjects),1);
                    corrTrials = cell(numel(subjects),1);
                    numEvent = cell(numel(subjects),1);
                    trigNum = cell(numel(subjects),1);
                    
                    for iTask = 1:size(task,1)
                        subData = get(task(iTask),'OutPutArguments');
                        subData = subData{1};
                        % combine subject specific results of preprocessing
                        subErp{iTask,1} = subData.subErp;
                        subErpEqual{iTask,1} = subData.subErpEqual;
                        subTrialInd{iTask,1} = subData.subTrialInd;
                        corrTrials{iTask,1} = subData.corrTrials;
                        numEvent{iTask,1} = subData.numEvent;
                        rejEpoch(iTask,:) = subData.rejEpoch;
                        trialNum{iTask,1} = subData.trialNum;                 
                        rejLog(iTask) = subData.rejLog;
                    end
                    % retrieve analysis paraemters from processed data set
                    analysis = subData.analysis;
                    % retrieve path parameters from processed data set
                    paths = subData.paths;
                    % set raw and result Dir
                    paths.rawDirAll = paths.local.rawDir;
                    paths.resDirAll = paths.local.resDir;
                    % retrieve trigger parameters from processed data set
                    trig = subData.trig;
                    % retrieve filter parameters from processed data set
                    filtPar = subData.filtPar;
                    dur.Start = get(job,'StartTime');
                    dur.End = get(job,'FinishTime');
                    save_erp(subErp,subErpEqual,subTrialInd,corrTrials,numEvent,rejEpoch,trialNum,subjects,paths,trig,dur,analysis,filtPar,rejLog,job);
                else
                    disp(':: There are unfinished tasks for the current job');
                    task
                end
              case 'plot'
                sched = [];
                job = [];
                % set raw and result Dir
                paths.rawDirAll = paths.local.rawDir;
                paths.resDirAll = paths.local.resDir;
              otherwise
                error([':: There is no option called ' parallelMode '.']);
            end
        end
    end % preprocessing
        
    if (~analysis.parallel && nargin < 3) || strcmpi(parallelMode,'plot')
        %% Calculation of difference
        % calculate difference waves
        [erpAll, restoredConf, chanlocs, trig] = calc_diff(subjects,paths,analysis,taskType);
        % Get Plot Parameters
        [plotPar] = config('Plot','task',taskType);
        % Select Channels to plot
        plot_erp(erpAll,chanlocs,plotPar,trig,analysis,paths,taskType,restoredConf);
    end

    
    