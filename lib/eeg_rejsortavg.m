% eeg_rejsortavg - Mark trials reducing SNR for rejection
%
% Usage:
%   >> EEG = eeg_rejsortavg(EEG, 'key1', value1, 'key2', value2, ...
%                                'keyn', valuen);
%
% Inputs:
%   EEG       - EEGLAB EEG structure
%
% Optional inputs:
%   'chans'   - vector channels to use for RMS and SP var {default all}
%   'pnts'    - vector sampling points to use for RMS and SP var {default
%               all}
%   'plot'    - (0|1) plot snr2 diagram {default 0}
%
% Output:
%   EEG       - EEGLAB EEG structure
%
% References:
%   Rahne, T., et al. (2008). Sorted averaging - application to auditory
%   event related responses. J Neurosci Meth, 172, 74-8
%
% Authors: Andreas Widmann and Alexandra Bendixen, University of Leipzig,
%   2010

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2010 Andreas Widmann and Alexandra Bendixen, University
% of Leipzig, widmann@uni-leipzig.de
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Id$

function [ EEG, sortIdx ] = eeg_rejsortavg(EEG, varargin)

Arg = struct(varargin{:});
if ~isfield(Arg, 'chans') || isempty(Arg.chans)
    Arg.chans = 1:size(EEG.data, 1);
end
if ~isfield(Arg, 'pnts') || isempty(Arg.pnts)
    Arg.pnts = 1:size(EEG.data, 2);
end

% Compute RMS
rmsArray = sqrt(mean(EEG.data(Arg.chans, Arg.pnts, :) .^ 2, 2)); % Time domain RMS
rmsArray = mean(rmsArray, 1); % Average RMS across channels

% Sort epochs by RMS
[foo sortIdx] = sort(rmsArray);
sortIdx = squeeze(sortIdx);
% sortIdx = 1:size(EEG.data, 3); % Debug only

% Compute SNR2
snr2 = zeros(1, size(EEG.data, 3));
sumSqArray = zeros(length(Arg.chans), length(Arg.pnts));
meanArray = EEG.data(Arg.chans, Arg.pnts, sortIdx(1));
for iEpoch = 2:size(EEG.data, 3)
    [sumSqArray meanArray] = rsumsq(EEG.data(Arg.chans, Arg.pnts, sortIdx(iEpoch)), sumSqArray, meanArray, iEpoch - 1);
    snr2(iEpoch) = iEpoch / mean(mean(sumSqArray / iEpoch, 2), 1);
end

% Plot SNR2 on request
if isfield(Arg, 'plot') && Arg.plot
    plot(snr2)
end

% Peak SNR2
[peak peakIdx] = max(snr2);

% Mark trials reducing SNR for rejection
if ~isfield(EEG, 'reject') || ~isfield(EEG.reject, 'rejmanual') || isempty(EEG.reject.rejmanual)
    EEG.reject.rejmanual = zeros(1, size(EEG.data, 3));
end
EEG.reject.rejmanual(sortIdx(peakIdx + 1:end)) = 1;

end

function [sumSqArray, meanArray, n] = rsumsq(xArray, sumSqArray, meanArray, n)
% recsumsq - Recursive computation of sum of squares
%
% References:
%   http://groups.google.ca/group/sci.math.stat/msg/26e6eb3ebaaa6558

    n = n + 1;
    devArray = xArray - meanArray;
    meanArray = meanArray + devArray / n;

    sumSqArray = sumSqArray + devArray .* (xArray - meanArray);

end