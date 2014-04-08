% Draw significant intervals in ERP plot.
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

function drawSig(sigInt,plotPar,analysis,timeRes)

allColor = {[1 0.8 0.8; 0.8 0.8 1];
            [1 0.3 0.3; 0.3 0.3 1]};

for iDraw = 1:2
    
    drawInt = zeros(1,size(sigInt.r,2));
    drawInt(2,:) = sigInt.r;

    % define color depending on loop iteration (light versions of the
    % colors in first iteration for non corrected data and dark versions
    % off the colors for the second iteration for multiple co mparison
    % corrected data)
    if iDraw == 1
        drawInt(1,sigInt.raw) = 1;
        color = allColor{1};
    else
        drawInt(1,sigInt.fdr) = 1;
        color = allColor{2};
    end
    

    % draw significant intervalls
    switch lower(plotPar.statTest)
      case 'anova'
        for iInt = 1:size(drawInt,2)
            if drawInt(1,iInt)
                intStart = iInt * timeRes + analysis.erpWin(1);
                % draw current intervall
                rectangle('Position',[intStart plotPar.yScale(2)-plotPar.yCoef*0.1 timeRes plotPar.yCoef*0.1],...
                          'FaceColor',color(1,:),'EdgeColor','none');
            end
        end
      case 'trendtest'
        for iInt = 1:size(drawInt,2)
            if drawInt(1,iInt)
                if drawInt(2,iInt) > 0
                    col = color(1,:);
                else
                    col = color(2,:);
                end
                intStart = iInt * timeRes + analysis.erpWin(1);
                % draw current intervall
                rectangle('Position',[intStart plotPar.yScale(2)-plotPar.yCoef*0.1 timeRes plotPar.yCoef*0.1],...
                          'FaceColor',col,'EdgeColor','none');
            end
        end
    end
end