% Check if all requested raw Data is available and prepare path
% specification 
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

function paths = checkRawData(paths,subjects,taskType)

% redefine paths to data according to task type
paths.rawDirAll = paths.rawDir;
paths.resDirAll = paths.resDir;
% get subjects with available raw data
paths.subDirs = dir([paths.rawDirAll paths.rawSubFolderPrefix '*']);

% get all subjects (all available raw data)
paths.allSubjects = zeros(1,size(dir([paths.rawDirAll paths.rawSubFolderPrefix '*']),1));

%% Get Numbers of all available subjects (available raw data)
for iSubj = 1:size(dir([paths.rawDirAll paths.rawSubFolderPrefix '*']),1)
    % for each raw data folder parse the subject number out of the folder name give
    % all the subject folders have the same subject specifier (e.g. vp) as given in 
    % the config (config.m)
    ind(iSubj).end = numel(paths.subDirs(iSubj).name);
    ind(iSubj).length = ind(iSubj).end - numel(paths.rawSubFolderPrefix);
    ind(iSubj).start = numel(paths.subDirs(iSubj).name) - (ind(iSubj).length - 1);

% find all subject numbers according to available folders
paths.allSubjects(1, iSubj) = str2num(paths.subDirs(iSubj).name(ind(iSubj).start:ind(iSubj).end));  
end

% take all subjects when no subjects are specified
if nargin < 2
    subjects = paths.allSubjects;
end

%% check if all requested subject data exists
if ~all(ismember(subjects,paths.allSubjects))
    missSubj = find(~ismember(subjects,paths.allSubjects));
    disp(' ');
    for iMsg = 1:numel(missSubj)
        disp([':: Did not find raw data for subject ' ...
              num2str(subjects(missSubj(iMsg))) '!']);
    end
    disp(' ');
    error(':: Some raw data is missing');
end
