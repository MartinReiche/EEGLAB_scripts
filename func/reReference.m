% Perform rereferencing to given reref channel in config.m
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

function erpAll = reReference(erpAll,chanlocs,analysis)
%% find new rererence channel labels and channel indices
% initialize array for rereference channel indices
    rerefInd = zeros(1,size(analysis.rerefChan,2));
    for iRefChan = 1:size(analysis.rerefChan,2)
        foundChan = 0;
        for iChan = 1:numel(chanlocs)
            if strcmp(analysis.rerefChan{iRefChan},chanlocs(iChan).labels)
                rerefInd(iRefChan) = iChan;
                foundChan = 1;
            end
        end
        if ~foundChan
            error([':: Did not find channel ' analysis.rerefChan{iRefChan} ' for re-referencing']);
        end
    end

    % select reference channels
    refChan = erpAll(:,:,:,rerefInd);
    
    if size(rerefInd,2) > 1
        % average given rereference channels
        refChan = mean(refChan,4);
    end

    disp(':: Re-referencing channels');
    % get channel list
    chanList = 1:size(erpAll,4);
    % % subtract each channel except mastoids

    
    % erpAll(:,:,:,chanList(~ismember(chanList,rerefInd))) = bsxfun(@minus,erpAll(:,:,:,chanList(~ismember(chanList,rerefInd))),refChan);
    % subtract each channel
    erpAll(:,:,:,chanList) = bsxfun(@minus,erpAll(:,:,:,chanList),refChan);
    
    