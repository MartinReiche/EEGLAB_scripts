% Returns given number of distinguishable colors by variing the "Hue" value
% in the HSV color space and keeping constant the "Stauration" and "Value"
% parameters. HSV values are then converted to RGB space.
% 
% Usage:
%   >> [colorMatrix] = distinct_color(n)
%                                     
% INPUT
% n                      - number of requested distinguishable colors
%
% OUTPUT
% colorMatrix            - matrix of n lines representing the colors and 3
%                          columns for the RGB values (between 0 and 1)
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

function [colorMatrix] = distinct_color(n)

% initialize color matrix
colorMatrix = zeros(n,3);
% define hue values depending on n
colorMatrix(:,1) = [1/n:1/n:1]';
colorMatrix(colorMatrix == 1) = 0;
% set saturation
colorMatrix(:,2) = 1;
% set value
colorMatrix(:,3) = 1;
% convert to rgb
colorMatrix = hsv2rgb(colorMatrix);
