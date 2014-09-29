% Handle the baseline correction if the rejection mode is other than sorted
% averaging. For data with sorted averaging, the baseline correction is
% already performed in segmentation.m

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

function erpAll = baseline_corr(erpAll,analysis,restoredConf)

if analysis.rmBase
    disp(':: Performing baseline correction'); 
    % get ms range of baseline window relative to beginning of epoch
    baseMS(1) = analysis.baseWin(1) - analysis.erpWin(1);
    baseMS(2) = analysis.baseWin(2) - analysis.erpWin(1);
    % correct for 0
    baseMS(baseMS == 0) = 1;
    % get the point range of the baseline window
    pointrange = ceil(baseMS(1)*analysis.sampRate/1000):floor(baseMS(2)*analysis.sampRate/1000);

    % separately for all subjects and all channels
    for iSubj = 1:size(erpAll,1)
        for iChan = 1:size(erpAll,4)
            % get the mean value in the baseline window of the current subject on the
            % current channel and subtract it from the whole data range of the
            % current subject on the current channel
            tmpmean = mean(double(erpAll(iSubj,:,pointrange,iChan)),3);
            erpAll(iSubj,:,:,iChan) = erpAll(iSubj,:,:,iChan) - repmat(tmpmean, [1 1 size(erpAll,3) 1]);
        end
    end

elseif  ismember(restoredConf.analysis.rejmode,[1 2 3]) && ~analysis.rmBase
    disp(':: WARNING. Skipping baseline correction!');
end