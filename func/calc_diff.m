% Calculate diffference waves for all triggers and all Subjects.
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


function [erpAll, restoredConf, chanlocs, trig] = calc_diff(subjects,paths,analysis,taskType)

%% Load ERP file
% get ERP file list
erpFile = dir([paths.resDirAll 'ERP(*']);

if size(erpFile,1) > 1
    
    disp(' ');
    disp(':: Detected more than one ERP file');
    disp(' ');
    for iFile = 1:size(erpFile,1)
        disp(['   File ' num2str(iFile) ': ' erpFile(iFile).name]);
    end
    disp(' ');
    disp(':: Choose one file (number)');
    iFile = input('>> ','s');
    % List of possible file indices
    posFiles = [1:size(erpFile,1)];

    % check validity of user input
    askAgain = 1;
    while askAgain
        if ismember(str2num(iFile),posFiles)
            askAgain = 0;
        else
            disp([':: ' num2str(iFile) ' is no possible file number']);
            disp([':: Choose one file (number) or type b - browse for additional files: ']);
            iFile = input('>> ','s'); 
            askAgain = 1; 
        end
    end

    disp(' ');
    disp(':: Loadind ERP file');
    disp(' ');
    if strcmpi(iFile,'b')
        uiopen('load')
    else
        load([paths.resDirAll erpFile(str2num(iFile)).name])
    end
elseif size(erpFile,1) == 1

    disp(' ');
    disp(':: There is one ERP file. Choose this file or browse for others (y/b)');
    disp(' ');
    disp(['   File: ' erpFile(1).name]);
    disp(' ');
    answ = input('>> ','s');
    askAgain = 1;
    while askAgain
        if any(strcmpi(answ,{'y','b'}))
            askAgain = 0;
        else
            disp([':: ' answ ' is no possible possible Option']);
            disp([':: Type y to load the file or  b - browse for additional files: ']);
            disp(' ');
            disp(['   File: ' erpFile(1).name]);
            disp(' ');
            answ = input('>> ','s'); 
            askAgain = 1; 
        end
    end
    
    if strcmpi(answ,'y')
        erpFile = dir([paths.resDirAll 'ERP(*']);
        disp(' ');
        disp(':: Loadind ERP file');
        disp(' ');
        load([paths.resDirAll erpFile(1).name])
    elseif strcmpi(answ,'b')
        uiopen('load')
        disp(' ');
        disp(':: Loadind ERP file');
        disp(' ');
    end

else
    disp('Did not find a ERP file, please select.');
    uiopen('load')
    disp(' ');
    disp(':: Loadind ERP file');
    disp(' ');
end

% retrieve erpAll from erp structure
if analysis.equalErp
    erpAll = erp.erpAllEqual; 
else
    erpAll = erp.erpAll;
end

% get chanlocs from erp structure
chanlocs = erp.chanlocs;
% restore triggers from timepoint of analysis
trig = erp.trig(end);

if analysis.savedErpWin 
    % restore config paramers from analysis
    restoredConf.analysis = erp.analysis(end);
    restoredConf.trig = erp.trig(end);
    restoredConf.paths = erp.paths(end);
    restoredConf.filtPar = erp.filtPar(end);
    % restore erpWindow
    analysis.erpWin = erp.analysis(end).erpWin;
else
    restoredConf = [];
end
clear erp;
%% Re-evaluate trigger colors
% get current trigger parameters
currTrigPars = config('triggers','task',taskType);


% build trigger and label arrays for given task
for iTrig = 1:size(trig.triggers,1)
    labelFound = 0;
    for iCurrTrig = 1:size(currTrigPars.triggers,1)
        if strcmpi(trig.triggers{iTrig,2},currTrigPars.triggers{iCurrTrig,2})
            labelMatch = iCurrTrig;
            trig.color(iTrig,:) = currTrigPars.color(iCurrTrig,:); 
        end
    end
end

%% determine indices of current trigger matrix for difference waves
trig.diffInd = {};
trig.diffCol = [];
iCount = 0;
iMin = [];
iSub = [];

for iDiff = 1:size(currTrigPars.diffWaves,1)
    for iTrig = 1:size(trig.trigLabels,1)
        switch trig.trigLabels{iTrig}
          case currTrigPars.diffWaves{iDiff,1}
            % find Index of minuend
            iMin = iTrig;
          case currTrigPars.diffWaves{iDiff,2}
            % find Index of subtrahend
            iSub = iTrig;
        end
    end
    if ~isempty(iMin) && ~isempty(iSub)
        % if Subtrahend and minuend was found for current difference wave
        % raise count variable
        iCount = iCount + 1; 
        % save index of minuend and subtrahend
        trig.diffInd{iCount,1} = [iMin iSub];
        % save the corresponding difference wave label
        trig.diffInd{iCount,2} = currTrigPars.diffWaves{iDiff,3};
        % clear minuend and subtrahend index variables
        iMin =  [];
        iSub = [];
        % get the color index of the current difference wave and store it in trigger
        % struct
        trig.diffCol(iCount,1:2) = cell2mat(currTrigPars.diffWaves(iDiff, 4:5));
    end
end

%% Compute difference waves

% Select current Subjects
erpAll = erpAll(subjects(:),:,:,:);
% rereference electrodes
if analysis.reref
    erpAll = reReference(erpAll,chanlocs,analysis);
end
% perform baseline correction
erpAll = baseline_corr(erpAll,analysis);

% Calculate difference waves
for iDiff = 1:size(trig.diffInd,1)
    erpAll(:,end+1,:,:) = erpAll(:,trig.diffInd{iDiff,1}(1),:,:)-erpAll(:,trig.diffInd{iDiff,1}(2),:,:);
    % Add index of current condition in erpAll structure in trig.diffInd for
    % later access
    trig.diffInd{iDiff,3} = size(erpAll,2);
end

% check whether new difference waves from diffference waves should be formed
trig.diffInd2 = {};
newCol = [ ];
iCount = 0;
iMin = [];
iSub = [];


for iDiff = 1:size(currTrigPars.diffWaves,1)

    for iTrig = 1:size(trig.diffInd,1)

        switch 1
            
          case strcmp(trig.diffInd{iTrig,2},currTrigPars.diffWaves{iDiff,1})
            % find Index of minuend
            iMin = trig.diffInd{iTrig,3};
            
          case strcmp(trig.diffInd{iTrig,2},currTrigPars.diffWaves{iDiff,2})
            % find Index of subtrahend
            iSub = trig.diffInd{iTrig,3};
            
        end
    end

    if ~isempty(iMin) && ~isempty(iSub)
        % if Subtrahend and minuend was found for current difference wave
        % raise count variable
        iCount = iCount + 1; 
        % save index of minuend and subtrahend
        trig.diffInd2{iCount,1} = [iMin iSub];
        % save the corresponding difference wave label
        trig.diffInd2{iCount,2} = currTrigPars.diffWaves{iDiff,3};
        % newCol
        % clear minuend and subtrahend index variables
        iMin =  [];
        iSub = [];
        % get the color index of the current difference wave and store it in trigger
        % struct
        trig.diffCol(iDiff,1:2) = cell2mat(currTrigPars.diffWaves(iDiff, 4:5));
    end
end

% Calculate difference waves
for iDiff = 1:size(trig.diffInd2,1)
    erpAll(:,end+1,:,:) = erpAll(:,trig.diffInd2{iDiff,1}(1),:,:)-erpAll(:,trig.diffInd2{iDiff,1}(2),:,:);
    % Add index of current condition in erpAll structure in trig.diffInd2 for
    % later access
    trig.diffInd2{iDiff,3} = size(erpAll,2);
end

% Concatenate new and old trig.diffInd
trig.diffInd = [trig.diffInd; trig.diffInd2];

