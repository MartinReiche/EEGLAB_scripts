% Takes vector of P-values and a given q-value and returns significance
% threshold based on false discovery rate
% 
% Author: unknown

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

function [pID,pN] = fdr(p,q)
% FORMAT [pID,pN] = fdr(p,q)
% 
% p   - vector of p-values
% q   - False Discovery Rate level
%
% pID - p-value threshold based on independence or positive dependence
% pN  - Nonparametric p-value threshold

% p = p(fininte(p));  % Toss NaN's
% modified AB because failed to find fininte 
if any(isnan(p))
    error('remove NaN entries from the p vector');
end
p = sort(p(:));
V = length(p);
I = (1:V)';

cVID = 1;
cVN = sum(1./(1:V));

pID = p(max(find(p<=I/V*q/cVID)));
pN = p(max(find(p<=I/V*q/cVN)));

