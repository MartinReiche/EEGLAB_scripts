% Part a given EEG structure for files specified in paths.partFile at
% sampling point specified in paths.partFile
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

function [EEG,paths] = part_file(EEG,partFile,paths)

for iPart = 1:size(partFile{1,4},1)
    ALLEEG(iPart) = pop_select( EEG,'point',partFile{1,4}(iPart,:));
end


if paths.partInd == 0
    % if this is the first part of the file
    % set the part index to the next file part (2 by default)
    paths.partInd = 2;
    % Assign first part of the file to EEG structure
    % EEG = ALLEEG(1);
    EEG = ALLEEG(1);
elseif paths.partInd == size(ALLEEG,2)
    % if this is the last part od the file
    % Assign first part of the file to EEG structure
    EEG = ALLEEG(size(ALLEEG,2));
    % terminate partition by setting paths.partInd to 0
    paths.partInd = 0;
else
    % Assign current part of the File to EEG structure
    EEG = ALLEEG(paths.partInd);
    % and increase the index by 1
    paths.partInd = paths.partInd + 1;
end

% Check for newly included 'boundary' events and remove them
endEvent = size(EEG.event,2);
iEvent = 1;
while iEvent <= endEvent
    % check all events sequentially
    if strcmpi(EEG.event(iEvent).type,'boundary')
        % if the current event is 'boundary',rm it
        disp(' ');
        disp([':: Removing boundary event for event index: ' num2str(iEvent)]);
        EEG.event(iEvent) = [];
        endEvent = endEvent - 1;
    else
        iEvent = iEvent + 1;
    end
end

disp(':: Reconvert events to numeric');
for iEvent = 1:size(EEG.event,2)
EEG.event(iEvent).type = str2num(EEG.event(iEvent).type);
end

