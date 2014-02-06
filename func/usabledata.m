% This functions asks the user which ERP file to read and calculates the
% amount of usable data for the subjects in the given file. A table is
% displayed with the amount of usable data per subject and subjects
% exceeding a given threshold are marked. A vector of subject numbers which
% are still within the acceptable bounds of usable data is returned.
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

function subjects = usabledata(paths,threshold)

%% Load ERP file
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

%% Calculate usable data
for iSubj=1:size(erp.erpAll,1)
    for iTrig = 1:size(erp.trig.triggers,1)
        % Get Proportion of usable Data for current subject and trigger
        corrProp(iSubj,iTrig) = erp.corrTrials{iSubj}(iTrig)/erp.eventCount{iSubj}(iTrig);
        if isnan(corrProp(iSubj,iTrig))
            corrProp(iSubj,iTrig) = 0; 
        end
    end
end
disp('------------------------------------------------------------------');
disp('                     USABLE DATA PER SUBJECT');
disp('------------------------------------------------------------------');
disp(['Rejection threshold: ' num2str(threshold, '%2.1f') ' %']);
disp(' ');
subjects = [];
for iSubj=1:size(erp.erpAll,1)

    currProp = round(mean(corrProp(iSubj,:))*1000)/10;
    if currProp < threshold
        disp([':: Subject ' num2str(iSubj, '%0.2d') '     ' num2str(currProp,'%2.1f') ' %      <- REJECT']);
    else
        disp([':: Subject ' num2str(iSubj, '%0.2d') '     ' num2str(currProp,'%2.1f') ' %']);
        subjects = [subjects iSubj];
    end
end
disp('------------------------------------------------------------------')
disp([':: Overall        ' num2str(round(mean(mean(corrProp,2))*1000)/10,'%2.1f') ' % (based on condition means)']);
disp('------------------------------------------------------------------')
end
