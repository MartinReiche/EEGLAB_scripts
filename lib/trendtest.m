% Perform within-subjects linear contrast analysis (linear trend test). Takes
% Data in matrix form where lines represent the subjects with the data of
% severeal factors represented by columns. Returns P-value for given data.
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

function p = trendtest(data)
    
% get number of levels
    nIV = size(data,2);
    % get number of samples
    n = size(data,1);
    
    % get linear coefficient weights
    coef = 1:nIV;
    coef = coef - mean(coef);

    % calculate sum of squares for contrast
    SSc = (n*(sum(coef.*mean(data)))^2)/sum(coef.^2);
    % calculate residual sum of squares for contrast
    SSce = (sum(sum((data .* repmat(coef,size(data,1),1)),2).^2) - ...
            (n*(sum(coef.*mean(data)))^2)) / sum(coef.^2);
    % calculate mean square of contrast residuals
    MSce = SSce/(n-1);
    % F ratio
    F = SSc/MSce;
    % obtain p-value from F distribution
    p = 1-fcdf(F,1,n-1);
    


