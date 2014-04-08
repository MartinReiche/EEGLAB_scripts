% Calculate plotting parameters (numbers of columns or rows depending on
% different input parameters) for example number of subjects or number of
% conditions
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

%% Get Display Arrangement parameters
function dispPars = getDisp(method,varargin)

% set defaults
plotLegend = 0;

%% evaluate input parameters
for iArg = 1:2:length(varargin)
   switch lower(varargin{iArg})
     case 'parameters'
       plotPar = varargin{iArg+1};
     case 'structure'
       nPlot = varargin{iArg+1};
     case 'legend'
       plotLegend = varargin{iArg+1};
     otherwise
       error([':: Wrong input argument: ' varargin{iArg}]);
   end
end

% validate input arguments
switch 1
  case ~exist('plotPar','var')
    error(':: Not enough input arguments, ''parameters'' is missing')
  case ~exist('nPlot','var')
    error(':: Not enough input arguments, ''structure'' is missing')
end

% switch method
switch lower(method)
  case 'electrode array'
    %% Display Arrangement of pre-spcecified electrode array
    if plotPar.autoDim
        % get number of plots
        numPlot = numel(nPlot);
        % add one plot slot for the legend
        if plotLegend
            numPlot = numPlot + 1;
        end
        % Determine grid arrangement depending on number of subjects
        dispPars.nCol = ceil(sqrt(numPlot));
        dispPars.nRow = ceil(numPlot/dispPars.nCol);
        dispPars.legend = dispPars.nRow * dispPars.nCol;
        dispPars.pos = [1:numel(nPlot)];
    else
        % take prespecified parameters from configfile
        dispPars.nRow = plotPar.plotDim(1);
        dispPars.nCol = plotPar.plotDim(2);
        dispPars.legend = plotPar.plotDim(3);
        dispPars.pos = plotPar.plotChannelPos;
    end
  case 'statistics'
    %% Display ERPs with Histograms (amplitude) at specific channel
    % get number of plots
    numPlot = numel(nPlot);
    % Determine grid arrangement depending on number of subjects
    dispPars.nCol = numPlot;
    dispPars.nRow = 2;
    dispPars.curvePos = [1:dispPars.nCol];
    dispPars.histPos = dispPars.curvePos + dispPars.nCol;
  otherwise
    error([':: There is no method called ' method]);
end






