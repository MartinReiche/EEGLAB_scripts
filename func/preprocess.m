% Perform preprocessing steps per subjects and return subject specific ERP
% data, and rejection related statistics
%
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

% Please define all parameters in config.m
function out = preprocess(taskType,iSubj,analysis,filtPar,trig,paths)

    tStart = tic; % start measuring time for Analysis

    %% run preprocessing routine
    % check availabel raw data and configure raw data paths
    paths = checkRawData(paths,iSubj,taskType);

    % prepare matrix for counting of available events for each subjects
    % (lines represent blocks, colums 1-20 are triggers 101-512)
    eventCount = zeros(size(trig.triggers,1),analysis.nBlocks);
    % prepare result folder for subject
    paths = prepSubDir(analysis,paths,iSubj);

    % get condition order
    condOrder = preporder(iSubj);
    counter = 0;
    
    for iFile = 1:analysis.nBlocks
        counter = counter + 1;
        % load raw data for current file and stimulation parameters
        [EEG, analysis, paths] = loadRawData(paths, iSubj,taskType,iFile,analysis,counter);
        % check parameters (sampling rate, triggers etc)
        EEG = checkFileBasic(EEG,iSubj,iFile,taskType,analysis,trig,paths,condOrder);
        % Retriggering and systematically exclude events from analysis
        EEG = change_trig(EEG,analysis,trig,iFile,condOrder,taskType); 
        % bipolarize eye channels
        EEG = bipolarize(EEG,analysis);
        % adjust electrode positions
        EEG = chanLoc(EEG,paths);
        % Run Channel interpolation
        EEG = interpChan(EEG,analysis,iSubj,iFile);
        % filter the block depending on rejection method:
        switch analysis.rejmode
            % No eye correction
          case {0,1,3,4}
            % filter Data
            filtPar.eye = 0; % no eye correction
            EEG = fir_filter(EEG,analysis,filtPar);
            
           % safe block data after segemtation and baseline correction
            eventCount = segmentation(EEG,trig,analysis,paths,eventCount, ...
                                      iFile,iSubj,condOrder);
            % eye correction
          case {2,5}
            % filter Data
            filtPar.eye = 1; % before eye correction
            EEG = fir_filter(EEG,analysis,filtPar);
            % save EEG block
            pop_saveset(EEG,[paths.resFileSubSpec num2str(iSubj,'%0.2d') ...
                             paths.resFileBlockSpec num2str(iFile, '%0.2d') '.set'], ...
                        paths.resDir);
          otherwise
            error([':: Invalid option for rejection method: ' num2str(analysis.rejmode)]);
        end
        clear EEG;
    end
    
    switch analysis.rejmode
        % no eye correction
      case {0,1,3,4}
        % merge files of same trigger for current subject and save
        merge_data(iSubj,paths,trig,eventCount,condOrder);
        % eye correction
      case {2,5}
        % merge all files for current subject, perform eye correction, apply post
        % filter and save
        merge_all(iSubj,paths,filtPar,analysis);
        % safe block data after trigger resegemtation and baseline correction
        segmentation([],trig,analysis,paths,eventCount, ...
                     iFile,iSubj,condOrder);
    end
    
    % artifact rejection and averaging
    [out.subErp,out.subErpEqual,...
     out.subTrialInd,out.corrTrials,...
     out.numEvent,out.rejEpoch,...
     out.trialNum,out.rejLog] = eeg_rejection(iSubj,paths,trig,analysis);
    
    % get timing 
    out.time = toc(tStart);
    
    % add config data to outputs for later retrieval
    out.analysis = analysis;
    out.filtPar = filtPar;
    out.trig = trig;
    out.paths = paths;

    % if this is the last subject of the current Job and the result dir of
    % the current job is empty, remove the result dir of the current job
    resFolders = dir([paths.resDirAll 'vp*']);
    if max(analysis.subjects) == iSubj && analysis.parallel && size(resFolders,1) == 0
        rmdir(paths.resDirAll);
    end
    
    % clear unused variables
    clear analysis filtPar trig paths eventCount condOrder;

