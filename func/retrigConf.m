% Takes current file and condition order and returns the retrigger Cell
% depending on the condition of the current file. Each trigger specified in
% the left column of the cell will be replaced with the trigger in the
% right column of the corresponding line of the retrigger cell.
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

function reTrig = retrigConf(iFile,condOrder)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!! Review change_trig.m before conducting a new analysis
% !!! To correct for new omission trigger range
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% retrigger events (first of pair/ second of pair)
reTrig = {
% NEW
% change certain proximal[+] distal[+]
    [41],          ((100 * condOrder(iFile)) + 41);  
% change certain proximal[+] distal[-]
    [42],          ((100 * condOrder(iFile)) + 42);  
% repetition possible proximal[+] distal[+]
    [43],          ((100 * condOrder(iFile)) + 43);  
% repetition possible proximal[+] distal[-]
    [44],          ((100 * condOrder(iFile)) + 44);  
% repetition possible proximal[-] distal[+]
    [45],          ((100 * condOrder(iFile)) + 45);  
% repetition possible proximal[-] distal[-]
    [46],          ((100 * condOrder(iFile)) + 46);  
% TENTATIVE
    [7],           ((100 * condOrder(iFile)) + 7);  
         };
end