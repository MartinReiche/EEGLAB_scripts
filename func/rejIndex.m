% Get indices for rejected trials per subject and trigger from former
% rejection procedure
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

function rejEpochs = rejIndex(paths)

rejEpochFile = dir([paths.resDirAll 'ERP(*']);

if size(rejEpochFile,1) > 1
    
    disp(' ');
    disp(':: Detected, more than one rejection file');
    disp(' ');
    for iFile = 1:size(rejEpochFile,1)
        disp(['   File ' num2str(iFile) ': ' rejEpochFile(iFile).name]);
    end
    disp(' ');
    disp(':: Choose one file (number)');
    iFile = input('>> ');
    % List of possible file indices
    posFiles = [1:size(rejEpochFile,1)];

    % check validity of user input
    askAgain = 1;
    while askAgain
        if ismember(iFile,posFiles)
            askAgain = 0;
        else
            disp([':: ' num2str(iFile) ' is no possible file number']);
            disp([':: Choose one file (number): ']);
            iFile = input('>> '); 
            askAgain = 1; 
        end
    end

    disp(' ');
    disp(':: Loadind rejection file');
    disp(' ');
    load([paths.resDirAll rejEpochFile(iFile).name])
elseif size(rejEpochFile,1) == 1
    rejEpochFile = dir([paths.resDirAll '*rejectedEpochs*']);
    disp(' ');
    disp(':: Loadind rejection file');
    disp(' ');
    load([paths.resDirAll rejEpochFile(1).name])
else         
    error(':: Rejection index file doesn''t exist');
end

rejEpochs = erp.rejEpochs;
clear erp;
