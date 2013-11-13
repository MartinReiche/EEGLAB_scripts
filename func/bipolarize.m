function EEG = bipolarize(EEG,analysis)
% Bipolarize upper and lower eye channel to vertical EOG and lateral eye
% channels to vertical EOG
%
% Author: Martin Reiche, martin.reiche@uni-oldnburg.de    
% based on a script by S.Andersen (2005)

%% find eye channels
disp(':: Replacing: HEOG = LO1-LO2, VEOG = SO1-IO1');
% initialize channel numbers
eyeChanNum = [];
for iChan = 1:size(EEG.chanlocs,2)
    for iEyeChan = 1:size(analysis.eyeChan,2)
        if strcmp(EEG.chanlocs(iChan).labels,analysis.eyeChan(:,iEyeChan))
            eyeChanNum(iEyeChan)=iChan;
        end
    end
end
%% Check if all eye channels were found and if perform bipolarization
if size(eyeChanNum,2)==size(analysis.eyeChan,2)
    EEG.data(eyeChanNum(1),:,:)=EEG.data(eyeChanNum(1),:,:)-EEG.data(eyeChanNum(2),:,:);
    EEG.chanlocs(eyeChanNum(1)).labels='HEOG';
    EEG.data(eyeChanNum(3),:,:)=EEG.data(eyeChanNum(3),:,:)-EEG.data(eyeChanNum(4),:,:);
    EEG.chanlocs(eyeChanNum(3)).labels='VEOG';
    EEG = pop_select(EEG,'nochannel',[eyeChanNum(2) eyeChanNum(4)]);
else
    error(':: Failed to find correct channels for subtraction')
end    
