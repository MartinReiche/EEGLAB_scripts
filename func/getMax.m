% Get the maximal y Scale values for one figure
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

function maxVal = getMax(method,erpAll,analysis,plotPar,plotConds,labels,iChan,channels2plot)

% evaluate method
switch lower(method)
  case 'electrode array'
    currInd = [];    

    % get curve indices
    for iCurve = 1:size(plotConds,2)
        % go through each Curve of the current figure
        foundLabel = 0;
        for iIndex = 1:size(labels,1)
            % check current curve label against each label of the
            % curves stored in erpAll to get the index of the current
            % curve label 
            if strcmp(plotConds(iCurve),labels(iIndex))
                foundLabel = 1;
                currInd = [currInd iIndex];
            end        
        end    
    end

    % get maxima for each curve on each channel
    for iChan = 1:size(channels2plot,2)
        for iWave = 1:numel(currInd)
            erpMax(iChan,iWave) = max(squeeze(mean(erpAll(:,currInd(iWave),:,channels2plot(iChan)))));
            erpMin(iChan,iWave) = min(squeeze(mean(erpAll(:,currInd(iWave),:,channels2plot(iChan)))));
            for iWin = 1:size(plotPar.compWin,1)
                statWin = [round(((plotPar.compWin(iWin,1)+abs(plotPar.xScale(1)))*analysis.sampRate)/1000) ...
                           round(((plotPar.compWin(iWin,2)+abs(plotPar.xScale(1)))*analysis.sampRate)/1000)];
                erpMean(iWin,iWave,iChan) = mean(mean(squeeze(erpAll(:,currInd(iWave),statWin,channels2plot(iChan))),2));
                erpErr(iWin,iWave,iChan) = std(mean(squeeze(erpAll(:,currInd(iWave),statWin,channels2plot(iChan))),2))/sqrt(size(erpAll,1));
            end
        end
    end

    % get maximal and minimal values (amplitude mean in time window plus SEM)
    meanMaxValues = erpMean + erpErr;
    meanMinValues = erpMean - erpErr;
    
    % assign output
    maxVal.erpMax = ceil(max(max(erpMax)));
    maxVal.erpMin = floor(min(min(erpMin)));
    maxVal.meanMin = min(min(min(squeeze(meanMinValues))));
    maxVal.meanMax = max(max(max(squeeze(meanMaxValues))));

  case 'statistics'
    % initialize variables
    currInd = [];
    
    for iPlot = 1:size(plotConds,2)
        % get the data to plot 
        for iCurve = 1:size(plotConds{iPlot},2) 
            % go through each Curve of the current figure
            foundLabel = 0;
            for iIndex = 1:size(labels,1)
                % check current curve label against each label of the
                % curves stored in erpAll to get the index of the current
                % curve label 
                if strcmp(plotConds{iPlot}(iCurve),labels(iIndex))
                    foundLabel = 1;
                    currInd = [currInd iIndex];
                end        
            end
        end
    end
    
    % get data
    for iWave = 1:numel(currInd)
        erpMax(iWave) = ceil(max(squeeze(mean(erpAll(:,currInd(iWave),:,channels2plot(iChan)),1))));
        erpMin(iWave) = floor(min(squeeze(mean(erpAll(:,currInd(iWave),:,channels2plot(iChan)),1))));
        for iWin = 1:size(plotPar.compWin,1)
            statWin = [round(((plotPar.compWin(iWin,1)+abs(plotPar.xScale(1)))*analysis.sampRate)/1000) ...
                       round(((plotPar.compWin(iWin,2)+abs(plotPar.xScale(1)))*analysis.sampRate)/1000)];
            erpMean(iWin,iWave) = mean(mean(squeeze(erpAll(:,currInd(iWave),statWin,channels2plot(iChan))),2));
            erpErr(iWin,iWave) = std(mean(squeeze(erpAll(:,currInd(iWave),statWin,channels2plot(iChan))),2))/sqrt(size(erpAll,1));
        end
    end 

    % get maximal and minimal values (amplitude mean in time window plus SEM)
    meanMaxValues = [];
    for iWave = 1:numel(erpMean)
        meanMaxValues = [meanMaxValues erpMean(iWave) - erpErr(iWave)];
        meanMaxValues = [meanMaxValues erpMean(iWave) + erpErr(iWave)];
    end

    % get minimal and maximal value for current figure
    maxVal.meanMin = min(meanMaxValues);
    maxVal.meanMax = max(meanMaxValues);
    maxVal.erpMin = min(erpMin);
    maxVal.erpMax = max(erpMax);
  
  otherwise
    error([':: Invalid Option: ' method]);
end