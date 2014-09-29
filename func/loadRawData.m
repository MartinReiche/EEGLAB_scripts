function [EEG, analysis, paths] = loadRawData(paths,nSubj,nSession,iFile,analysis,counter)
% script to load raw EEG data file using EEGLAB. Currently loads BIOSIG and
% Brainvision Recorder raw files. Needs to be adjusted for other raw data
% types.
%
% Copyright (c) 2013 Martin Reiche, Carl-von-Ossietzky-University Oldenburg
% Author: Martin Reiche, martin.reiche@uni-oldnburg.de
    
%% Parameters
% if there is another part of fileCat, decrease file increment to load the
% same raw file as last time
    
    for iPart = 1:size(paths.partFile,1)
        if (nSubj == paths.partFile{iPart,1}) & nSession == paths.partFile{iPart,2}
            rawFileOrder = paths.partFile{iPart,5};
        else
            rawFileOrder = 1:size(paths.allFiles,1);
    end
    rawFile = [paths.rawDir paths.allFiles{rawFileOrder(iFile)}];
    
    % last check if file exists
    currentFile = dir(rawFile);
    if isempty(currentFile)
        error([':: File ' rawFile ' dows not exist!']);
    end
    clear curentFile;
    
    %% get channel number of nose for current subject
    % if iFile == 1
    
    if counter == 1 && strcmpi(analysis.rawFormat,'biosemi')
        foundRef = 0;
        disp(' ');
        disp([':: Get number of reference channel for subject ' num2str(nSubj, ...
                                                          '%0.2d')]);
        disp(' ');
        EEG = pop_biosig(rawFile);
        for iChan = 1:numel(EEG.chanlocs)
            if strcmpi(EEG.chanlocs(iChan).labels,analysis.refChan)
                analysis.refChanNum = iChan;
                foundRef = 1;
            end
        end
        if foundRef
            disp(' ');
            disp([':: Found Reference (' analysis.refChan ') at channel ' ...
                  num2str(iChan) ' for subject ' num2str(nSubj,'%0.2d')]);
        else
            error([':: Did not find Reference (' analysis.refChan ') for subject ' ...
                   num2str(nSubj,'%0.2d')]);
        end
    end

    %% Load current File
    disp(' ');
    disp('############################################');
    disp(' ');
    disp([':: Loading raw data of Subject ' num2str(nSubj,'%0.2d') ' Block ' num2str(analysis.blocks(iFile),'%0.2d')]);
    disp(' ');
    disp('############################################');
    disp(' ');

    % evaluate raw data format and call respective loading routine
    switch analysis.rawFormat
      case 'biosemi'
        % biosemi system
        EEG = pop_biosig(rawFile,'ref',analysis.refChanNum);
      case 'brainvision'
        % brainvision analyzer files
        EEG = pop_loadbv(paths.rawDir,paths.allFiles{rawFileOrder(iFile)}, [], []);
    end

    % check if reference channel number is right
    if strcmpi(analysis.rawFormat,'biosemi')
        if ~strcmpi(EEG.chanlocs(1).ref,analysis.refChan) 
            input([':: Detected wrong reference channel. Expected: ' analysis.refChan ', got: ' ...
                   EEG.chanlocs(1).ref]);
        else
            % remove reference channel
            disp(' ');
            disp(':: Removing reference channel');
            disp(' ');
            EEG = pop_select(EEG,'nochannel',analysis.refChanNum);
        end
    end

    for iPart = 1:size(paths.partFile,1)
        if nSubj == paths.partFile{iPart,1} & nSession == paths.partFile{iPart,2} & rawFileOrder(iFile) == paths.partFile{iPart,3}
            % part file
            [EEG,paths] = part_file(EEG,paths.partFile(iPart,:),paths);
        end
    end

    % delete channel EXG8
    foundExg = 0;
    for iChan = 1:numel(EEG.chanlocs)
        if strcmpi(EEG.chanlocs(iChan).labels,analysis.exgChan)
            exgChanNum = iChan;
            foundExg = 1;
        end
    end
    if foundExg
        disp(' ');
        disp(':: Deleting EXG8 channel');
        disp(' ');
        EEG = pop_select(EEG,'nochannel',[exgChanNum]);
    elseif ~isempty(analysis.exgChan)
        error([':: Did not find channel: ' analysis.exgChan ' for subject ' ...
               num2str(nSubj,'%0.2d')]);
    end

    % % go through the chaninterp array
    % % load corresponding stimulation parameter file
    % load([...
    %     paths.behavDir...
    %     paths.behavSubjSpec num2str(nSubj,'%0.2d')...
    %     paths.behavBlockSpec num2str(iFile,'%0.2d')...
    %     paths.behavFileExt]);

end
