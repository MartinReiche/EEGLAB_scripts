% Perform running one-way ANOVA over each time point at given Channel and
% return significance index over the whole time range. Output vector
% [sigInt] contains indices of datapoints which reach significance
% according to fdr with diven q value (plotPar.alpha) in config.m
% 
% USAGE: sigInt = runningStat(erpAll,plotPar)
%
% Copyright (C) 2013 Martin Reiche, Carl-von-Ossietzky-University Oldenburg
% Author: Martin Reiche, martin.reiche@uni-oldnburg.de

% NEEDS:
%
% ANOVA Scripts by Trujillo-Ortiz, A., R. Hernandez-Walls and
%      R.A. Trujillo-Perez. (2004). RMAOV1:One-way repeated measures ANOVA. A
%      MATLAB file. [WWW document]. URL
%      http://www.mathworks.com/matlabcentral/
%      fileexchange/loadFile.do?objectId=5576

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

function sigInt = runningStat(erpAll,spData)

%% determine parameters
% perform the statistics only when the statistics time range is within the bounds of the epoch
if (spData.plotPar.runStatWin(1) < spData.analysis.erpWin(1)) || (spData.plotPar.runStatWin(2) > spData.analysis.erpWin(2))
    disp([':: Statistics time range exceeds bounds of epoch. Testing whole epoch range']);
    spData.plotPar.runStatWin(1) = spData.analysis.erpWin(1);
    spData.plotPar.runStatWin(2) = spData.analysis.erpWin(2);
elseif (spData.plotPar.runStatWin(1) > spData.analysis.erpWin(1)) || (spData.plotPar.runStatWin(2) < spData.analysis.erpWin(2))
    disp(' ');
    warning(':: NOT TESTING WHOLE EPOCH RANGE!');
    disp(' ');
end

% get time resolution
timeRes = 1000/spData.analysis.sampRate;
% get relevant range for testing (in ms from beginning
statRangeMS = [abs(spData.analysis.erpWin(1) - spData.plotPar.runStatWin(1))...
               abs(spData.analysis.erpWin(1) - spData.plotPar.runStatWin(2))];

% convert ms timerange to datapoints
statRange = [ceil(statRangeMS(1)/timeRes) floor(statRangeMS(2)/timeRes)];
% convert zeros to one
statRange(statRange == 0) = 1;
% get relevant data
statData = erpAll(:,spData.currInd,statRange(1):statRange(2),spData.channelIndex);
% initialize array to store p values
pVals = zeros(1,size(statData,3));    
% Go through all the timepoints
for iPoint = 1:size(statData,3)
    % perfrom one-way RMANOVA with the factor CONDITION (available waves per plot) 
    pVals(iPoint) = OneWayrmAoV(statData(:,:,iPoint));
    
end    
pID = fdr(pVals,spData.plotPar.alpha);
if ~isempty(pID)
    sigInt = find(pVals < pID);
else
    sigInt = [];
end
