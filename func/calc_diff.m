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
% Select current Subjects
erpAll = erpAll(subjects(:),:,:,:);
% rereference electrodes
if analysis.reref && ~analysis.gfp
    erpAll = reReference(erpAll,chanlocs,analysis);
end
% perform baseline correction
% check rejection mode and force baseline correction on

if ismember(restoredConf.analysis.rejmode,[4 5]) && ~analysis.rmBase 
    % automatically enable baseline correction
    disp([':: Forcing baseline correction for ' analysis.rejLabel{restoredConf.analysis.rejmode + 1}]); 
    analysis.rmBase = 1;
end
erpAll = baseline_corr(erpAll,analysis,restoredConf);



% build trigger and label arrays for given task
for iTrig = 1:size(trig.triggers,1)
    labelFound = 0;
    for iCurrTrig = 1:size(currTrigPars.triggers,1)
        if strcmpi(trig.triggers{iTrig,2},currTrigPars.triggers{iCurrTrig,2})
            trig.color(iTrig,:) = currTrigPars.color(iCurrTrig,:); 
        end
    end
end
trig.diffWaves
if ~isempty(trig.diffWaves)
    %% determine indices of current trigger matrix for difference waves
    % get matrix of difference waves to calculate 
    toCalc = currTrigPars.diffLabels;
    disp(':: Get difference wave index');
    while ~isempty(toCalc)
        origTrig = trig.trigLabels;
        [trig,toCalc] = getDiffInd(currTrigPars,trig,toCalc);
        if isequal(origTrig,trig.trigLabels)
            for iCurve = 1:size(toCalc,1)
                disp([':: Missing curve label for difference wave: ' num2str(toCalc{iCurve})]);
            end
            error('Did not find curve labels for calculation of difference waves.') 
        end
    end

    % Calculate difference waves
    disp(':: Calculate difference waves');
    for iDiff = 1:size(trig.diffInd,1)
        erpAll(:,end+1,:,:) = erpAll(:,trig.diffInd{iDiff,1}(1),:,:)-erpAll(:,trig.diffInd{iDiff,1}(2),:,:);
        % Add index of current condition in erpAll structure in trig.diffInd for
        % later access
        trig.diffInd{iDiff,3} = size(erpAll,2);
    end


    % Calculate difference waves
    for iDiff = 1:size(trig.diffInd,1)
        erpAll(:,end+1,:,:) = erpAll(:,trig.diffInd{iDiff,1}(1),:,:)-erpAll(:,trig.diffInd{iDiff,1}(2),:,:);
        % Add index of current condition in erpAll structure in trig.diffInd for
        % later access
        trig.diffInd{iDiff,3} = size(erpAll,2);
    end
    trig.color = [trig.color; trig.diffCol];

end


function [trig,toCalc] = getDiffInd(currTrigPars,trig,toCalc)

% initialize fields in trigger structure for difference wave calculation
    if ~isfield(trig,'diffInd')
        trig.diffInd = {};
        trig.diffCol = [];
        trig.iCount = 0;
    end
    % initialize minuend and subtrahend variables
    iMin = [];
    iSub = [];
    calculatedDiff = {};

    for iDiff = 1:size(currTrigPars.diffLabels,1)
        if ismember(currTrigPars.diffLabels{iDiff},toCalc)
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
                trig.iCount = trig.iCount + 1; 
                % save index of minuend and subtrahend
                trig.diffInd{trig.iCount,1} = [iMin iSub];
                % save the corresponding difference wave label
                trig.diffInd{trig.iCount,2} = currTrigPars.diffWaves{iDiff,3};
                % clear minuend and subtrahend index variables
                iMin =  [];
                iSub = [];
                % get the color index of the current difference wave and store it in trigger
                % struct
                trig.diffCol(trig.iCount,1:2) = cell2mat(currTrigPars.diffWaves(iDiff, 4:5));
                % delete current difference label in toCalc 
                index = find(strcmp(toCalc, currTrigPars.diffLabels{iDiff}));
                toCalc(index) = [];
                % get name of calculated difference wave
                calculatedDiff = [calculatedDiff; currTrigPars.diffLabels{iDiff}];
            else
                % clear minuend and subtrahend index variables
                iMin =  [];
                iSub = [];
            end
        end
    end

    % Add Calculated Diff names to trigger labels
    trig.trigLabels = [trig.trigLabels; calculatedDiff];

