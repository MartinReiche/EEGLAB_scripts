% Prepare raw and result dirs for each subject and check whether there
% already is any data which mighrt interfere with the current analysis procedure

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

function paths = prepSubDir(analysis,paths,iSubj)

foundDir = 0;
for iDir = 1:size(paths.subDirs,1)
   if strcmp([paths.rawSubFolderPrefix num2str(iSubj, '%0.2d')],paths.subDirs(iDir).name);
       subDirInd = iDir;
       foundDir = 1;
   end       
end

if ~foundDir
    error([':: Did not find raw data folder for subject ' num2str(iSubj)]);
end

% change results and raw paths according to subject
paths.rawDir = [paths.rawDirAll paths.subDirs(subDirInd).name '/'];
paths.resDir = [paths.resDirAll paths.subDirs(subDirInd).name '/'];

% construct file names fpr current subject
paths.allFiles = cell(numel(analysis.blocks),1);
for iFile = 1:numel(analysis.blocks)
    paths.allFiles{iFile} = [paths.rawSubFilePrefix num2str(iSubj,'%0.2d')...
                             paths.rawFileBlockSpec num2str(analysis.blocks(iFile),'%0.2d')... 
                             paths.rawFileExt];
end

% check for old epoch files extracted from separate blocks and delete them
if size(dir([paths.resDir ...
             paths.resFileSubSpec '*' ... 
             paths.resFileBlockSpec '*' ...
             paths.resFileTrigSpec '*.set']),1)
    % if some old epoch files already existed, delete them
    disp(' ');
    disp([':: Detected old epoch files for subject ' ...
          num2str(iSubj,'%0.2d') ', deleting them']);
    delete([paths.resDir '*' paths.resFileBlockSpec '*']);
end

% check whether the subjects results folder already exists and if not
% create it
if ~exist(paths.resDir)
    disp([':: Creating folder for subject ' num2str(iSubj, '%0.2d')]);
    mkdir(paths.resDir);
end