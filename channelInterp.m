% determine interpolation parameters for given channels in given subjects
% for passive and active task

function out = channelInterp(taskType)

% Subject, channel name, block
% PASSIVE TASK
chanInterp1 = {
% EXAMPLE: 

% interpolate channel number 62 of subject 3 for trigger files containing
% triggers 501, 503 and 511    
% 3, 62, [501 503 511]; 

% interpolate channel 'P4' of subject 6 in raw file number 16
% 6, 'P4', 16;

% interpolate channel 'CPz' of subject 9 in raw file number 1 to 18
% 9, 'CPz', [1:18];        
              };
% ACTIVE TASK
chanInterp2 = {
% interpolation configuration for 2nd task
              };

% define outputs
switch taskType
  case 1
    out = chanInterp1;
  case 2
    out = chanInterp2;
  otherwise
    error([':: Task type ' num2str(taskType) ' is not available']); 
end