% Merge all block data in one file (to prepare for eye correction)
% and apply low-pass filter
%
% Copyright (c) 2013 Martin Reiche, Carl-von-Ossietzky-University Oldenburg
% Author: Martin Reiche, martin.reiche@uni-oldnburg.de

function EEG = merge_all(iSubj,paths,filtPar,analysis)

% get list of all files
allFiles = dir([paths.resDir...
                paths.resFileSubSpec num2str(iSubj,'%0.2d')...
                paths.resFileBlockSpec '*.set']);

% load all preprocessed block files
for iFile = 1:size(allFiles,1)
    ALLEEG(iFile) = pop_loadset(allFiles(iFile).name,paths.resDir);
end
EEG = pop_mergeset(ALLEEG,[1:size(allFiles,1)],0);

% Find EOG Channels
eogVector = [];
for iChannel = 1:size(EEG.data,1) %determine relevant channel numbers
    if strcmp(EEG.chanlocs(iChannel).labels,'HEOG') || strcmp(EEG.chanlocs(iChannel).labels,'VEOG')
        eogVector = [eogVector iChannel];
    end
end
if numel(eogVector) ~= 2
    disp(':: Failed to find eye channels');
    error(':: Try again');
end

%eye movement correction
disp(' ');
disp(':: Performing eye correction.');
EEG = eeg_emcp(EEG,'eeg',setdiff(1:size(EEG.data,1),eogVector),'eog',eogVector,'avgsub',false);


% filter Data
filtPar.eye = 2; % after eye correction
EEG = fir_filter(EEG,analysis,filtPar);

% save Data
EEG = pop_saveset(EEG,[paths.resFileSubSpec num2str(iSubj,'%0.2d')...
                    'all.set'],paths.resDir);

% Delete old files .set & .fdt
disp(':: Deleting old block Files');
delete([paths.resDir... 
        paths.resFileSubSpec  num2str(iSubj,'%0.2d')...
        paths.resFileBlockSpec '*.set']);
delete([paths.resDir... 
        paths.resFileSubSpec  num2str(iSubj,'%0.2d')...
        paths.resFileBlockSpec '*.fdt']);
% Clear EEG data structures to prevent memory overflow
clear ALLCOM ALLEEG CURRENTSET EEG LASTCOM;
end