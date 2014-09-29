% Configuration file for EEG Analysis, should hold all the configurable
% parameters which than get called by the specific functions
% Usage:
%   >> [Argout 1, Argout 2, Argout n] = config(Method, 'key1', value1, 'key2', ...
%                                                   value2, 'keyn', valuen);
% INPUT
% Method              - Define which Parameters to return:
%
% OUTPUT
% Variable List of Paramters depending on Input Method
% Ex: Input Method 'Path' returns all 
% Triggers are defined in ./triggerlabels.m
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

function varargout = config(Method,varargin)
%% ANALYSIS PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run preprocessing
    analysis.preprocess = 1;
    % perform parallel preprocessing 
    analysis.parallel = 0;
    % Cluster configuration (enter name of zour cluster configuration here
    % [as under Parallel > Select Configuration])
    analysis.core = 'Cluster_configuration_name';
    % filter the data (logic switch for all filter routines)
    analysis.filterFlag = 1;
    % raw data format (options: 'biosemi', 'brainvision')
    analysis.rawFormat = 'brainvision';
    % change triggers according triggerlabels('retrig',trig,taskType)
    analysis.changeTrig = 0;
    % Eye channels for bipolarization, specify in following order
    % {'LO1', 'LO2', 'SO1', 'IO1'}
    analysis.eyeChan = {'LO1', 'LO2', 'SO1', 'IO1'};
    % Rejection mode (0: no rejection, 1 delta rejection, 2 delta + eye
    % correction, 3 reject events specified in file, 4 sorted averaging,
    % 5 sorted averaging  + eye correction)
    analysis.rejmode = 2;
    % perform baseline correction (this option is only available for
    % rejection modes other than sorted averaging, with sorted averaging
    % baseline correction will always be performed
    analysis.rmBase = 1;
    % enable/disable rereferencing
    analysis.reref = 1;
    % new reference channel (if several electrodes are given here, than the
    % new reference will be the average of the given electrodes)
    analysis.rerefChan = {'M1' 'M2'};
    % erp window
    analysis.erpWin = [-100 200];
    % baseline window 
    analysis.baseWin = [-100 0];
    % block numbers for passive stimulation
    analysis.blocksPassive = 1:16;
    % block numbers for active stimulation
    analysis.blocksActive = 17:20;
    % Indices of Events to exclude 
    analysis.eventExcl = [1 2];
    % Trigger of events which were systematically excluded
    analysis.exclTrig = '98';
    % Time window for correct resonse after omission (in ms)
    analysis.respWin = [50 1050];
    % Original response Trigger (for exclusion of trigger around them in
    % predefined range [analysis.respEx])
    analysis.respTrig = 'S 17';
    % deviant / omission triggers (for exclusion of trigger around them in
    % predefined range [analysis.omEx])
    analysis.omissionRange = {'S 13' 'S 14' 'S 23' 'S 24' 'S 33' 'S 34' 'S 43'...
                        'S 44'};
    % Sampling Rate of raw files (in Hz)
    analysis.sampRate = 500;
    % Range of events to exclude around response (in ms)
    analysis.respEx = [-310 410];   
    % Range of events to exclude around omission (in ms)
    analysis.omEx = [0 610]; 
    % use same amount of trials for standards and deviants
    analysis.equalErp = 0;
    % plot ERPs
    analysis.plotERPflag = 1;
    % Plot with Stats
    analysis.statsFlag = 1;
    % Plot Topographies
    analysis.topoFlag = 0; 
    % Electrode Name for Reference
    analysis.refChan = 'Nose';
    % Electrode Name of EXG Channel (will be excluded)
    analysis.exgChan = [];
    % Delta Criterion (in microVolt)
    analysis.sortthresh = 100;
    % maximal proportion of events rejected on one electrode to allow before
    % printing Electrode in rejFile
    analysis.chanMaxRej = 0.1;
    % Reject Epochs without activity
    analysis.rejFlatepochs = 1;
    % minimum voltage change per trial
    analysis.flatthresh = 0.2; 
    % delete subject folder after rejection and averaging
    analysis.clearFolders = 1;
    % restore erpWin from saved config parameters of erp file
    analysis.savedErpWin = 1;
    % Rejection method labels
    analysis.rejLabel = {'no rejection';'delta rejection';...
                        'delta rejection and eye-movement correction';...
                        'rejection based on predefined indices';...
                        'sorted averaging';'sorted averaging and eye-movement correction'};
    analysis.rejFileLabel = {'no_rej';'delta_rej';'delta_eye';'file_rej';'sorted_avr';'sorted_eye'};

    %% CHANNEL INTERPOLATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Subject, channel name, block
    % PASSIVE TASK
    chanInterp1 = channelInterp(1);

    % ACTIVE TASK
    chanInterp2 = channelInterp(2);

    % channels to exclude from analysis
    analysis.excludeElecs = {
        {'01';'none'}; %subject 01
        {'02';'none'}; %subject 02
        {'03';'none'}; %subject 03
        {'04';'none'}; %subject 04
        {'05';'none'}; %subject 05
        {'06';'none'}; %subject 06
        {'07';'none'}; %subject 07
        {'08';'none'}; %subject 08
        {'09';'none'}; %subject 09
        {'10';'none'}; %subject 10
        {'11';'none'}; %subject 11
        {'12';'none'}; %subject 12
        {'13';'none'}; %subject 13
        {'14';'none'}; %subject 14
        {'15';'none'}; %subject 15
        {'16';'none'}; %subject 16
        {'17';'none'}; %subject 17
        {'18';'none'}; %subject 18
        {'19';'none'}; %subject 19
        {'20';'none'}; %subject 20
                   };
           
    %% PATH CONFIGURATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Path to raw eeg data
    paths.local.rawDir = '/local/path/to/your/raw/files/';
    paths.remote.rawDir = '/path/to/raw/files/on/cluster/';
    % Path to results dir
    paths.local.resDir = '/local/path/to/save/result/files/';
    paths.remote.resDir = '/path/cluster/to/save/result/files';
    % task folder names (gets appenden to arw file destination path)
    paths.taskLabel = {'passive/' 'active/'};
    % Path to behavioral results
    paths.local.behavDir = '/local/path/to/stimulation/results/';
    paths.remote.behavDir = '/path/to/stimulation/results/on/cluster/';
    % Path to topographie figures 
    paths.topoDir = 'topographies/';
    % Path to stim functions
    paths.local.stimFuncDir = '/local/path/to/stimultion/functions/';
    paths.remote.stimFuncDir = '/path/to/stimulation/functions/on/cluster/';
    % Path to function dir and lib dir (contents get uploaded when analysis is
    % carried out on cluster)
    paths.funcDir = '/local/path/to/eeg/functions/';
    % Path to Analysis lib Dir
    paths.local.libDir = '/local/path/to/eeg/libraries/';
    paths.remote.libDir = '/path/to/eeg/libraries/on/cluster/';
    % EEGLAB Dir
    paths.local.eeglabDir = '/local/path/to/eeglab/';
    paths.remote.eeglabDir = '/path/to/eeglab/on/cluster/';
    % electrode stup file (specify only the file which is in the /lib folder
    % without absolut path [with lib as the root path])
    paths.elecSetup = 'elec_96ch.elp';
    % result file extension
    paths.resFileExt = '.set';
    % result subject folder prefix (names of result data subject folders)
    paths.resSubFolderPrefix = 'vp';
    % result file subject specifier
    paths.resFileSubSpec = 'Subj';
    % result file block specifiern
    paths.resFileBlockSpec = 'block';
    % result file trigger specifier
    paths.resFileTrigSpec = 'tr';
    % raw subject folder prefix (names of raw data subject folders)
    paths.rawSubFolderPrefix = 'vp';
    % raw data file type
    paths.rawFileExt = 'vhdr';
    % behavioral file subject specifier 
    paths.behavSubjSpec = 'Subj';
    % behavioral file block specifier 
    paths.behavBlockSpec = 'Block';
    % behavioral file extension
    paths.behavFileExt = '.mat';
    % concatenate file (1st: sub, 2nd: tasktype, 3rd: file number, sampling
    % point of border, raw file loading sequence)
    % EXAMPLE:
    % paths.partFile = {17, 2, 11, [1 93764; 93800 189276], [1:11 11 12:17]};
    paths.partFile = {[],[],[],[],[]};

    %% FILTER CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % settings WITHOUT/AFTER eye correction
    % filter designed with firfilt 1.5.5
    % 512 Hz sampling rate
    % 0.1 - 30 Hz passband 
    % 0.001 Max passband daviation/ripple (= -60 dB)
    % 1.0 transition bandwidth
    % window type: kaiser
    % Kaiser beta: 5.65326 
    % Filter order: 1856
    % filter display name
    filtPar.name = '48 - 52 Hz notch filter';
    % filter type
    filtPar.fType = 'bandstop';
    % filter pass band
    filtPar.pass = [48 52];
    % filter order
    filtPar.fOrder = 1856;
    % window type
    filtPar.wType = 'kaiser';
    % Kaiser beta
    filtPar.kaiserBeta = 5.65326;
    % % % % % % % % % % % % % % % % % % % % % % %

    % settings BEFORE eye correction
    % filter designed with firfilt 1.5.5
    % 512 Hz sampling rate
    % 0.1 - 100 Hz passband
    % 0.001 passband deviation/ripple (= -60 dB)
    % 1.0 transition bandwidth
    % window type: kaiser
    % Kaiser beta: 5.65326
    % Filter order: 1856

    % enable the pre filter
    filtPar.pre.enable = 1;
    % filter display name
    filtPar.pre.name = '0.1 - 100 Hz band pass filter';
    % filter type
    filtPar.pre.fType = 'bandpass';
    % filter pass band
    filtPar.pre.pass = [0.1 100];
    % filter order
    filtPar.pre.fOrder = 1856;
    % window type
    filtPar.pre.wType = 'kaiser';
    % Kaiser beta
    filtPar.pre.kaiserBeta = 5.65326;
    % % % % % % % % % % % % % % % % % % % % % % %
    
    % settings AFTER  eye correction
    % filter designed with firfilt 1.5.5
    % 512 Hz sampling rate
    % 30 Hz lowpass
    % 0.001 passband deviation/ripple (= -60 dB)
    % 1.0 transition bandwidth
    % window type: kaiser
    % Kaiser beta: 5.65326
    % Filter order: 1856

    % enable the post filter
    filtPar.post.enable = 1;
    % filter display name
    filtPar.post.name = '48 Hz low pass filter';
    % filter type
    filtPar.post.fType = 'lowpass';
    % filter pass band
    filtPar.post.pass = [48];
    % filter order
    filtPar.post.fOrder = 1856;
    % window type
    filtPar.post.wType = 'kaiser';
    % Kaiser beta
    filtPar.post.kaiserBeta = 5.65326;


    %% GET TRIGGER DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Configure triggers in triggerlabels.m
    
    %% PLOTTING CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Components with stats
    plotPar.comps = {'Win1'};
    % Time range of components
    plotPar.compWin = [-20 0];
    % draw the baseline interval 
    plotPar.drawBaseLine = 1;
    % baseline window
    plotPar.baseWin = analysis.baseWin;
    % Line width of ERP Graphs
    plotPar.lineWidth = 1.5;
    % Time scaling
    plotPar.xScale = analysis.erpWin;
    % amplitude scaling (auto scale if empty)
    plotPar.yScale = [];
    % amplitude scaling for bar diagram (auto scale if empty)
    plotPar.yScaleBar = [];
    % run point by point RMANOVA in each plot
    plotPar.runningStat = 1;
    % define test type for running statistics ('anova' or 'trendtest')
    plotPar.statTest = 'trendtest';
    % define alpha (q) level for fdr of running anova 
    plotPar.alpha = 0.05;
    % define time window for running statistics
    % plotPar.runStatWin = analysis.erpWin;
    plotPar.runStatWin = analysis.erpWin;
    % horizontal axis coefficient
    plotPar.xCoef = 100;
    % vertical axis coefficient
    plotPar.yCoef = 1;
    % use grid for electrode array plotting
    plotPar.grid = 0;
    % display grid for called single ERP plots
    plotPar.singleGrid = 1;
    % Channels to plot without stats
    plotPar.plotChannels = {'HEOG';
                    'VEOG';
                    'E16';
                    'E08';
                    'E09';
                    'E07';
                    'E02';
                    'E03';
                    'E41';
                    'E01';
                    'E38';
                    'E22';
                    'E05';
                    'E18'};
    % automatic plotting dimensions, if this is 0 the
    % dimensionsspecification below will be used
    plotPar.autoDim = 0;
    % automatically choose color for graphs in one subplot
    plotPar.autoColor = 0;
    % Channel Positions for Plotting
    plotPar.plotChannelPos = [1 3:15];
    % Dimensions for Plotting
    plotPar.plotDim = [5 3 2];
    % plotting of stimuli
    plotPar.drawStim = 0;
    % latencies and position on yScale (voltage) for stimuli drawings
    plotPar.stim = [-900 0.5; -750 -1.5; -600 -1.5; -450 -0.5; -300 -0.5; -150 1]; 
    % tone duration for drawn stimuli (in ms)
    plotPar.stimDur = 50;
    % automatic scaling of clicked single ERP (0: no - boundarys of whole
    % figure will be used, 1: yes - boundaries will be determined depending
    % on single ERP data)
    plotPar.singleScaleAuto = 0;
    % channels to exclude from topographies
    plotPar.noTopoChan = {};
    % Conditions to plot, each line of the cell array represents one
    % figure, in each figure all the curves specified in one cell of
    % plotConds is plotted at a specified set of electrodes, the last
    % column in each cell specifies in which task the curves of the current
    % cell occur
    plotConds{1} = {'first-tone-1','second-tone-1','first-omission-1','second-omission-1',[1 2]};

    % Color Setting
    % 1 = same color per condition [across condition comparison]
    % 2 = same color per type [within condition comparison]
    plotPar.plotCondsCol = 1;
    
    plotPar.plotConds = plotConds;
    %% Plot with Stats and Histogram
    % Overhead of voltage scaling (gets added to max values per figure)
    plotPar.yOverhead = 0.5;
    % Channels to plot with stats (per line)
    plotPar.plotChannelsStat = {'E02'};
    % plot topographies in color or grayscale
    plotPar.colorFlag = 1;
    % Figures to plot. Labels of curves within one Cell go in one
    % plot. Lines of the cell array represent separate plots in one figure
    % and columns represesent separate figures
    
    plotPar.plotCondsStat{1,1} = {'first-tone-5','first-tone-4','first-tone-3','first-tone-2','first-tone-1',[1 2]};
    plotPar.plotCondsStat{2,1} = {'second-tone-5','second-tone-4','second-tone-3','second-tone-2','second-tone-1',[1 2]};
    plotPar.plotCondsStat{3,1} = {'tone-diff-5','tone-diff-4','tone-diff-3','tone-diff-2','tone-diff-1',[1 2]};

    % Color Setting
    % 1 = same color per type [within condition comparison]
    % 2 = same color per condition [across condition comparison]
    plotPar.plotCondsStatCol = [2];
    
    % End OF PARAMETER ADJUSTMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
% Evaluate input Method and return Parameters
switch lower(Method)
  case 'parameters'
    % Define Output
    if strcmpi(varargin{1},'task')
        taskType = varargin{2};
    else
        error([':: There is no input argument called ' varargin{1}]);
    end
    
    if ~strcmpi(varargin{1},'task')
        error([':: Wrong option: ' varargin{1} ' (task is required)']);
    elseif taskType == 1
        analysis.chanInterp = chanInterp1; 
        analysis.blocks = analysis.blocksPassive;
        analysis.nBlocks = numel(analysis.blocks);
    elseif taskType == 2
        analysis.chanInterp = chanInterp2; 
        analysis.blocks = analysis.blocksActive;
        analysis.nBlocks = numel(analysis.blocks);
    else
        error([':: Invalid value for Option: ' varargin{1}])
    end

    % define output arguments
    varargout{1} = analysis;
    
  case 'path'
    taskFlag = 0;
    analysisFlag = 0;
    filtFlag = 0;
    for iArg = 1:2:length(varargin)
        switch lower(varargin{iArg})
          case 'task'
            taskType = varargin{iArg + 1};
            taskFlag = 1;
          case 'analysis'
            analysis = varargin{iArg + 1};
            analysisFlag = 1;
          case 'filt'
            filtPar = varargin{iArg + 1};
            filtFlag = 1;
          otherwise
            error([':: Invalid Option ' varargin{iArg} ' for input Method ' Method]);
        end
    end
    if ~taskFlag | ~analysisFlag
       error([':: ''task'', ''filt'' and ''analysis'' is required for input Method ''' Method '''']) 
    end

    % determine paths according to cluster configuration (local or remote)
    if ismember(analysis.core,{'local.singleCore','local.dualCore', ...
                            'local.tripleCore','local.quadCore'}) || ~analysis.parallel
        paths.rawDir = paths.local.rawDir;
        paths.resDir = paths.local.resDir;
        paths.behavDir = paths.local.behavDir;
        paths.eeglabDir = paths.local.eeglabDir;
        paths.libDir = paths.local.libDir;
        paths.stimFuncDir = paths.local.stimFuncDir;
        % path to electrode setup file
        paths.elecSetup = [paths.libDir paths.elecSetup];
        paths.cluster = 'local';
    elseif ismember(analysis.core,{'HERO'})
        paths.rawDir = paths.remote.rawDir;
        paths.resDir = paths.remote.resDir;
        paths.behavDir = paths.remote.behavDir;
        paths.eeglabDir = paths.remote.eeglabDir;
        paths.libDir = paths.remote.libDir;
        paths.stimFuncDir = paths.remote.stimFuncDir;
        % path to electrode setup file
        paths.elecSetup = [paths.libDir paths.elecSetup];
        paths.cluster = 'remote';
    else
        error([':: Cluster configuration mismatch, could not configure ' ...
               'paths.']);
    end
    paths.chanlocs = paths.resDir;
    paths.rawDir = [paths.rawDir paths.taskLabel{taskType}];
    paths.behavDir = [paths.behavDir paths.taskLabel{taskType}];
    paths.local.topoDir = [paths.local.resDir paths.taskLabel{taskType} paths.topoDir];
    
    % add pre and post filter to rejFileLabel
    if ~filtPar.pre.enable
        analysis.rejFileLabel{3} = [analysis.rejFileLabel{3} '_noPreFilt'];
    end

    % define ERP and rejectedEpoch file names
    if ismember(analysis.rejmode,[0 1 3 4])
        if analysis.filterFlag
            if size(filtPar.pass,2) == 1
                filtName = [num2str(filtPar.pass) 'Hz_'];
            else
                filtName = [num2str(filtPar.pass(1)) '-' num2str(filtPar.pass(2)) 'Hz_'];
            end
        else
            filtName = 'NOFILT_';
        end
        paths.erpFileName = ['ERP(' num2str(analysis.erpWin(1)) 'to' num2str(analysis.erpWin(2)) 'ms)_'...
                            filtName analysis.rejFileLabel{analysis.rejmode+1} '.mat'];
        paths.rejFileName = ['rejectedEpochs(' num2str(analysis.erpWin(1)) 'to' num2str(analysis.erpWin(2)) ...
                            'ms)_' filtName analysis.rejFileLabel{analysis.rejmode+1} '.mat'];
    elseif ismember(analysis.rejmode,[2 5])

        if ~analysis.filterFlag | ~filtPar.post.enable
            filtName = 'NOFILT_';
        else
            if size(filtPar.post.pass,2) == 1
                filtName = [num2str(filtPar.post.pass) 'Hz_'];
            else
                filtName = [num2str(filtPar.post.pass(1)) '-' num2str(filtPar.post.pass(2)) 'Hz_'];
            end
        end

        paths.erpFileName = ['ERP(' num2str(analysis.erpWin(1)) 'to' num2str(analysis.erpWin(2)) 'ms)_'...
                            filtName analysis.rejFileLabel{analysis.rejmode+1} '.mat'];
        paths.rejFileName = ['rejectedEpochs(' num2str(analysis.erpWin(1)) 'to' num2str(analysis.erpWin(2)) ...
                            'ms)_' filtName analysis.rejFileLabel{analysis.rejmode+1} '.mat'];
    else
        error(':: Unknown rejection method');
    end

    % Define Output
    paths.partInd = 0;
    varargout = {paths};
    
  case 'filter'
    % Define Output
    varargout = {filtPar};
    
  case 'triggers'
    % evaluate input arguments
    for iArg = 1:2:numel(varargin)
        switch lower(varargin{iArg})
          case 'task'
            taskType = varargin{iArg+1};
          otherwise 
            error([':: There is no input argument called ' varargin{iArg}]);
        end
    end
    trig = [];
    % fetch trigger codes, trigger labels and diffwave configurations
    trig = triggerlabels('triggers',trig,taskType);
    % Initialize Trigger Matrix and Label array for given Task
    trig.newTriggers = {};
    trig.trigLabels = {};
    trig.color = [];
    if ~isempty(trig.diffWaves)
        trig.diffLabels = trig.diffWaves(:,3);
    else
        trig.diffLabels = [];
    end
    % build trigger and label arrays for given task
    for iTrig = 1:size(trig.triggers,1)
        if ismember(taskType,trig.triggers{iTrig,3})
            trig.newTriggers = {trig.newTriggers; trig.triggers{iTrig,1}};
            trig.trigLabels = [trig.trigLabels; trig.triggers{iTrig,2}];
            trig.color = [trig.color; trig.triggers{iTrig,4:5}];
        end
    end

    % % Throw Error because label not found
    % error(':: Problem finding labels, please Check spelling')    

    varargout = {trig};

  case 'plot'
    if strcmpi(varargin{1},'task')
        taskType = varargin{2};
    else
        error([':: There is no input argument called ' varargin{iArg}]);
    end
    
    if size(plotPar.comps,1) == size(plotPar.compWin,1)
        
        % Define Window markers for bar plots
        for iWin = 1:size(plotPar.comps,1)
            plotPar.winNames{iWin} = ['Window ' num2str(iWin)];
        end
        % Define Output Arguments with stats
        varargout = {plotPar};
    else
        error([':: plotPar.comps and plotPar.compWin do not have ' ...
               'the same size'])
    end

    % check for completeness of arguments
    if size(plotPar.plotCondsStat,2) ~= size(plotPar.plotCondsStatCol,2)
       error(':: ''plotCondsStat'' and ''plotCondsStatCol'' dimension Mismatch. Please adjust!.'); 
    end
    
    % Define Output Arguments without stats
    varargout = {plotPar};
    
  otherwise
    error([':: There is no input Method called ' Method]);
    
end

