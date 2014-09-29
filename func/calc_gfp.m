% Calculate the global field power of a given EEG time series on a given
% set of electrodes.
%
% USAGE: 
%
% 
% 
% 
% 
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

function gfp = calc_gfp(erpAll,chanlocs)

% get HEOG and VEOG position
eogHpos = strmatch('HEOG',{chanlocs.labels},'exact');
eogVpos = strmatch('VEOG',{chanlocs.labels},'exact');
if isempty(eogHpos) || isempty(eogVpos)
   error(':: Could not find HEOG and VEOG channels'); 
end

%get the data
erpTemp = erpAll(:,:,:,setdiff(1:size(erpAll,4),[eogHpos eogVpos])); %eliminate eye channels
nElecs = size(erpTemp,4);
nCond = size(erpTemp,2);
nSubj = size(erpTemp,1);
erpGFP = zeros(size(erpTemp,2),size(erpTemp,1),size(erpTemp,3),nElecs); %condition x subject x sampling point x channel
for iCond = 1:nCond
    erpGFP(iCond,:,:,:) = squeeze(erpTemp(:,iCond,:,:));
end

% %compute average reference
% for iCond = 1:nCond
%     for iSubj = 1:nSubj
%         erpTemp = squeeze(erpGFP(iCond,iSubj,:,:));
%         avgref = mean(erpTemp,2);
%         for iElec = 1:nElecs
%             erpGFP(iCond,iSubj,:,iElec) = erpTemp(:,iElec) - avgref;
%         end
%     end
% end

%compute GFP
gfp = zeros(nCond,nSubj,size(erpAll,3)); %condition x subject x sampling point
for iCond = 1:nCond
    for iSubj = 1:nSubj
        for iTime = 1:size(erpAll,3)
            erpTemp = squeeze(erpGFP(iCond,iSubj,iTime,:));
            gfp(iCond,iSubj,iTime) = sqrt(1/nElecs*sum(erpTemp(:).^2)); %or just use std(erpTemp) which gives the same result!
        end
    end
end

gfp = permute(gfp,[2 1 3]); 