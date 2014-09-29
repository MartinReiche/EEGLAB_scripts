% Set cannel locations to desired coordinates. Specify the coordinate file
% in config.m under paths.elecSetup
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

function EEG = chanLoc(EEG,paths)

% rename RM and LM to M2 and M1
for iChan = 1:numel(EEG.chanlocs)
    % go through all channels
    if strcmp(EEG.chanlocs(iChan).labels,'LM')
        EEG.chanlocs(iChan).labels = 'M1';
    end
    if strcmp(EEG.chanlocs(iChan).labels,'RM')
        EEG.chanlocs(iChan).labels = 'M2';
    end
end

% Set channel locations 
[paths.elecSetup]
EEG=pop_chanedit(EEG,'lookup',paths.elecSetup);
disp(' ');
end