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

function drawSig(sigInt,plotPar,analysis,timeRes,color)

% draw significant intervals
if ~isempty(sigInt)
    % find start and end of significant intervals
    foundStart = 0;
    intCount = 0;
    diffSig = [diff(sigInt) == 1 0];
    for iInt = 1:size(diffSig,2)
        if diffSig(1,iInt) && ~foundStart
            % increase the interval counter
            intCount = intCount + 1;
            % save start of current interval
            intBound(intCount,1) = sigInt(iInt);
            foundStart = 1;
        elseif ~diffSig(1,iInt) && foundStart
            % save end of current interval
            intBound(intCount,2) = sigInt(iInt);
            foundStart = 0;
        end
    end
    % draw significant intervalls
    for iInt = 1:size(intBound,1)

        % for each intervall
        % determine start end end of current intervall in ms relative to epoch start
        intStart = intBound(iInt,1) * timeRes + analysis.erpWin(1);
        intDur = (intBound(iInt,2) * timeRes +  analysis.erpWin(1))-intStart;
        % draw current intervall
        rectangle('Position',[intStart plotPar.yScale(2)-plotPar.yCoef*0.1 intDur plotPar.yCoef*0.1],...
                  'FaceColor',color,'EdgeColor','none');
    end
end

