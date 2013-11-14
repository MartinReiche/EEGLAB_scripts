% Batch analysis 
function batch_analysis(joblist)
taskType = {};
subjects = {};

% Seconds to wait bewteen Blocks
waitBlock = 5;

%% JOB 1 - Passive - Delta + EYE; PREFILTER; 45 Hz LowPass
jobIndex = 1;
% Cofiguration
subjects{jobIndex} = [1:20];
taskType{jobIndex} = 1;
% Initialize Configuration (DEFAULTS)
[analysis{jobIndex} filtPar{jobIndex} trig{jobIndex}] = call_default(taskType{jobIndex});
% CUSTOMIZATION
analysis{jobIndex}.rejmode = 2;
% FILTER
analysis{jobIndex}.filterFlag = 1;
filtPar{jobIndex}.pre.enable = 1;
filtPar{jobIndex}.pre.name = '0.1 - 100 Hz band pass filter';
filtPar{jobIndex}.pre.fType = 'bandpass';
filtPar{jobIndex}.pre.pass = [0.1 100];
filtPar{jobIndex}.pre.fOrder = 1856;
filtPar{jobIndex}.pre.wType = 'kaiser';
filtPar{jobIndex}.pre.kaiserBeta = 5.65326;
filtPar{jobIndex}.post.enable = 1;
filtPar{jobIndex}.post.name = '48 Hz low pass filter';
filtPar{jobIndex}.post.fType = 'lowpass';
filtPar{jobIndex}.post.pass = 48;
filtPar{jobIndex}.post.fOrder = 1856;
filtPar{jobIndex}.post.wType = 'kaiser';
filtPar{jobIndex}.post.kaiserBeta = 5.65326;
% ERP
analysis{jobIndex}.erpWin = [-300 0];
% CORE
analysis{jobIndex}.jobIndex = jobIndex;
analysis{jobIndex}.core = 'local.dualCore';

%% JOB 2 - Active - Delta + EYE; PREFILTER; 45 Hz LowPass
jobIndex = 2;
% Cofiguration
subjects{jobIndex} = [1:20];
taskType{jobIndex} = 2;
% Initialize Configuration (DEFAULTS)
[analysis{jobIndex} filtPar{jobIndex} trig{jobIndex}] = call_default(taskType{jobIndex});
% CUSTOMIZATION
analysis{jobIndex}.rejmode = 2;
% FILTER
analysis{jobIndex}.filterFlag = 1;
filtPar{jobIndex}.pre.enable = 1;
filtPar{jobIndex}.pre.name = '0.1 - 100 Hz band pass filter';
filtPar{jobIndex}.pre.fType = 'bandpass';
filtPar{jobIndex}.pre.pass = [0.1 100];
filtPar{jobIndex}.pre.fOrder = 1856;
filtPar{jobIndex}.pre.wType = 'kaiser';
filtPar{jobIndex}.pre.kaiserBeta = 5.65326;
filtPar{jobIndex}.post.enable = 1;
filtPar{jobIndex}.post.name = '48 Hz low pass filter';
filtPar{jobIndex}.post.fType = 'lowpass';
filtPar{jobIndex}.post.pass = 48;
filtPar{jobIndex}.post.fOrder = 1856;
filtPar{jobIndex}.post.wType = 'kaiser';
filtPar{jobIndex}.post.kaiserBeta = 5.65326;
% ERP
analysis{jobIndex}.erpWin = [-300 0];
% CORE
analysis{jobIndex}.jobIndex = jobIndex;
analysis{jobIndex}.core = 'local.dualCore';



% END OF JOB ASSIGNMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Run Batch analysis
for iJob = 1:numel(joblist)
    disp(' ');
    disp(' ');
    disp([' ### Starting Job ' num2str(joblist(iJob)) ' in ' num2str(waitBlock) ' seconds! ###']);
    disp(' ');
    disp(' ');
    pause(waitBlock);
    % Run current Job
    eeg_analysis(...
    taskType{joblist(iJob)},...
    subjects{joblist(iJob)},...
    analysis{joblist(iJob)},...
    filtPar{joblist(iJob)},...
    trig{joblist(iJob)},1);
end

end

function [analysis filtPar trig] = call_default(taskType)
% call default configuration structures for current job
analysis = config('parameters','task',taskType);
filtPar = config('filter');
trig = config('triggers','task',taskType);
end