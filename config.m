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
    analysis.parallel = 1;
    % Cluster configuration 'local...' ('.singleCore' '.dualCore' ...)
    % ex 'local.tripleCore'
    analysis.core = 'local.dualCore';
    % filter the data
    analysis.filterFlag = 1;
    % raw data format (options: 'biosemi', 'brainvision')
    analysis.rawFormat = 'brainvision';
    % change triggers according to func/retrigConf.m
    analysis.changeTrig = 0;
    % Eye channels for bipolarization, specify in following order
    % {'LO1', 'LO2', 'SO1', 'IO1'}
    analysis.eyeChan = {'E23', 'E17', 'E84', 'E29'};
    % Rejection mode (0: no rejection, 1 delta rejection, 2 delta + eye
    % correction, 3 reject events specified in file, 4 sorted averaging
    % [not yet implemented])
    analysis.rejmode = 2;
    % perform baseline correction (this option is only available for
    % rejection modes other than sorted averaging, with sorted averaging
    % baseline correction will always be performed
    analysis.rmBase = 0;
    % erp window
    analysis.erpWin = [-300 0];
    % baseline window 
    analysis.baseWin = [-300 -200];
    % Number of Blocks
    analysis.nBlocks = 15;
    % Indices of Events to exclude 
    analysis.eventExcl = [1 2];
    % Trigger of events which were systematically excluded
    analysis.exclTrig = '98';
    % Time window for correct resonse after omission (in ms)
    analysis.respWin = [50 1050];
    % Original response Trigger
    analysis.respTrig = 'S 17';
    % Sampling Rate of raw files (in Hz)
    analysis.sampRate = 500;
    % Range of events to exclude around response (in ms)
    analysis.respEx = [-310 410];   
    % Range of events to exclude around omission (in ms)
    analysis.omEx = [0 610]; % OLD [-160 610];   
    % use same amount of trials for standards and deviants
    analysis.equalErp = 0;
    % plot ERPs
    analysis.plotERPflag = 0;
    % Plot with Stats
    analysis.statsFlag = 1;
    % Plot Topographies
    analysis.topoFlag = 0; 
    % Electrode Name for Reference
    analysis.refChan = 'Nose';
    % Electrode Name of EXG Channel
    analysis.exgChan = [];
    % Delta Criterion
    analysis.sortthresh = 100;
    % maximal proportion of events rejected on one electrode to allow before
    % printing Electrode in rejFile
    analysis.chanMaxRej = 0.2;
    % Reject Epochs without activity
    analysis.rejFlatepochs = 1;
    % minimum voltage change per trial
    analysis.flatthresh = 0.2; 
    % delete subject folder after rejection and averaging
    analysis.clearFolders = 1;
    % restore erpWin from saved config parameters of erp file
    analysis.savedErpWin = 1;
    % Rejection method
    analysis.rejLabel = {'no rejection';'delta rejection';...
                        'delta rejection and eye correction';...
                        'rejection based on predefined indices';...
                        'sorted averaging'};
    analysis.rejFileLabel = {'no_rej';'delta_rej';'delta_eye';'file_rej';'sorted_avr'};

    %% CHANNEL INTERPOLATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Subject, channel name, block
    % PASSIVE TASK
    chanInterp1 = {
        % 3, 62, [501 503 511];
        % 4, 45, [301 303];      
        % 6, 'P4', 16;
        % 9, 'CPz', [1:18];        
        % 9, 'P4', [13:14];
        % 9, 'Oz', [13:14];
        % 9, 29, [102 112 401 402 403 411 412];
        % 12, 'P2', [17:18];
                  };
    % ACTIVE TASK
    % chanInterp2 = {
    %     1, 'F8', [1:18];
    %     3, 30, [311 312 321 412];
    %               };

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
    paths.rawDir = '/home/martin/documents/valax/raw_eeg/';
    % Path to results dir
    paths.resDir = '/home/martin/documents/valax/results/';
    % task folder names
    paths.taskLabel = {'passive/' 'active/'};
    % Path to behavioral results
    paths.behavDir = '/home/martin/Dropbox/PhD/valax/stim/results/';
    % Path to stim functions
    paths.topoDir = 'topographies/';
    paths.stimFuncDir = '/home/martin/Dropbox/PhD/valax/stim/func/';
    % Path to function dir and lib dir
    paths.funcDir = [pwd '/func/'];
    % Path to Analysis lib Dir
    paths.libDir = [pwd '/lib/'];
    % EEGLAB Dir
    paths.eeglabDir = '/home/martin/build/matlab11b/toolbox/eeglab12_0_2_5b/';
    % path to electrode setup file
    paths.elecSetup = [paths.libDir 'elec_96ch.elp'];
    % result file extension
    paths.resFileExt = '.set';
    % result subject folder prefix (names of result data subject folders)
    paths.resSubFolderPrefix = 'vp';
    % result file subject specifier
    paths.resFileSubSpec = 'Subj';
    % result file block specifier
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
    plotPar.compWin = [-30 0];
    % draw the baseline interval 
    plotPar.drawBaseLine = 1;
    % baseline window
    plotPar.baseWin = analysis.baseWin;
    % Line width of ERP Graphs
    plotPar.lineWidth = 1.5;
    % Time scaling
    plotPar.xScale = analysis.erpWin;
    % run point by point RMANOVA in each plot
    plotPar.runningStat = 1;
    % define time window for running statistics
    plotPar.runStatWin = [-300 0];
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
                    'F3';
                    'Fz';
                    'F4';
                    'FC1';
                    'FCz';
                    'FC2';
                    'C3';
                    'Cz';
                    'C4';
                    'M1';
                    'Pz';
                    'M2'};
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
    % Control of behavior in which figure the called subplot will be shown
    % (1 - open extra figure which gets overridden every time a new subplot
    % is called, 2 - open a new figure for each call)
    plotPar.singleDispMode = 1;
    % Conditions to plot, each line of the cell array represents one
    % figure, in each figure all the curves specified in one cell of
    % plotConds is plotted at a specified set of electrodes, the last
    % column in each cell specifies in which task the curves of the current
    % cell occur
    plotConds{1} = { 'Tone-RepEx-2','Tone-RepEx-3', 'Diff', [1 2]};
    plotConds{2} = { 'Tone-RepEx-4', [1]};

    % Color Setting
    % 1 = same color per condition [across condition comparison]
    % 2 = same color per type [within condition comparison]
    plotPar.plotCondsCol = 1;
    
    plotPar.plotConds = plotConds;
    %% Plot with Stats and Histogram
    % Overhead of voltage scaling (gets added to max values per figure)
    plotPar.yOverhead = 0.5;
    % Channels to plot with stats (per line)
    plotPar.plotChannelsStat = {'E01'};
    % plot topographies in color or grayscale
    plotPar.colorFlag = 1;
    % Figures to plot. Labels of curves within one Cell go in one
    % plot. Lines of the cell array represent separate plots in one figure
    % and columns represesent separate figures
    
    plotPar.plotCondsStat{1,1} = {'first-tone-1','first-tone-2','first-tone-3','first-tone-4','first-tone-5',[1 2]};
    plotPar.plotCondsStat{2,1} = {'second-tone-1','second-tone-2','second-tone-3','second-tone-4','second-tone-5',[1 2]};
    plotPar.plotCondsStat{3,1} = {'tone-diff-1','tone-diff-2','tone-diff-3','tone-diff-4','tone-diff-5',[1 2]};
    
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
    elseif taskType == 2
        analysis.chanInterp = chanInterp2; 
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
            error([':: Invalid Option ' varargin{iArg} ' for input Metho ' Method]);
        end
    end
    if ~taskFlag | ~analysisFlag
       error([':: ''task'', ''filt'' and ''analysis'' is required for input Method ''' Method '''']) 
    end
    
    paths.behavDir = [paths.behavDir paths.taskLabel{taskType}];
    paths.topoDir = [paths.resDir paths.taskLabel{taskType} paths.topoDir];
    
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
    elseif analysis.rejmode == 2 

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
    % Initialize label flag
    trig.plotLabels = 0;
    trig.diffLabels = 0;
    % evaluate input arguments
    for iArg = 1:2:numel(varargin)
        switch lower(varargin{iArg})
          case 'task'
            taskType = varargin{iArg+1};
          otherwise 
            error([':: There is no input argument called ' varargin{iArg}]);
        end
    end
    % fetch trigger codes, trigger labels and diffwave configurations
    trig = triggerlabels(trig,taskType);
    % Initialize Trigger Matrix and Label array for given Task
    trig.newTriggers = {};
    trig.trigLabels = {};
    trig.color = [];
    % build trigger and label arrays for given task
    for iTrig = 1:size(trig.triggers,1)
        if ismember(taskType,trig.triggers{iTrig,3})
            trig.newTriggers = [trig.newTriggers; trig.triggers{iTrig,1}];
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

