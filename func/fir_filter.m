% Filter EEG Data for given block using firfilt 1.5  by Andreas Widmann
% http://www.uni-leipzig.de/~biocog/content/widmann/eeglab-plugins/#firfilt

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


function [EEG] = fir_filter(EEG,analysis,filtPar)
    if analysis.filterFlag
        switch filtPar.eye
          case 0 % no eye correction
            % filter the data
            disp(' ' );disp([':: Apply ' filtPar.name]);disp(' ');
            [EEG,com,b] = pop_firws(EEG,'ftype',filtPar.fType,'fcutoff',filtPar.pass, ...
                                    'forder',filtPar.fOrder, 'wtype',filtPar.wType, ...
                                    'warg',filtPar.kaiserBeta);
            
            
          case 1 % before eye correction
            if filtPar.pre.enable
                disp(' ');disp([':: Apply ' filtPar.pre.name ' (PRE-FILTER)']);disp(' ');
                [EEG,com,b] = pop_firws(EEG,'ftype',filtPar.pre.fType,'fcutoff',filtPar.pre.pass, ...
                                        'forder',filtPar.pre.fOrder, 'wtype',filtPar.pre.wType, ...
                                        'warg',filtPar.pre.kaiserBeta);
            else
                disp(' ');
                disp('   !!! NO PRE FILTER IN USE !!! ');
                disp(' ');
            end
          case 2 % after eye correction
            if filtPar.post.enable
            disp(' ');disp([':: Apply ' filtPar.post.name ' (POST-FILTER)']);disp(' ');
            [EEG,com,b] = pop_firws(EEG,'ftype',filtPar.post.fType,'fcutoff',filtPar.post.pass, ...
                                    'forder',filtPar.post.fOrder, 'wtype',filtPar.post.wType, ...
                                    'warg',filtPar.post.kaiserBeta);
            else
                disp(' ');
                disp('   !!! NO POST FILTER IN USE !!! ');
                disp(' ');
            end
          otherwise
            error([':: Invalid Option ' num2str(analysis.rejmode) ' for filtPar.eyeRej']);
        end
        
    else
        disp(' ');
        disp('   !!! NO FILTER IN USE !!! ');
        disp(' ');
    end
end

