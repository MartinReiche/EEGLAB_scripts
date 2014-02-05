% Names of triggers and trigger configurations
function trig = triggerlabels(method,trig,taskType)
%% TRIGGER CONFIGURATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original Triggers
% 1 - first tone of pair (no preceding pair)
% 2 - second tone of pair (no preceding pair)
% 3 - first tone of pair (preceding pair)
% 4 - second tone of pair (preceding pair)
% 5 - tone of none pair (with preceding pair)
% 6 - tone of none pair (without preceding pair)
% 11 - omission @ first tone of pair (no preceding pair)
% 12 - omission @ second tone of pair (no preceding pair)
% 13 - omission @ first tone of pair (preceding pair)
% 14 - omission @ second tone of pair (preceding pair)
% 15 - omission @ tone of none pair (with preceding pair)
% 16 - omission @ tone of none pair (without preceding pair)

%% Control Triggers
% accepted latency threshold for stimuli (in ms)
trig.rejThresh.stim = 1; 
% accepted latency threshold for responses (in ms)
trig.rejThresh.resp = 5; 
% accepted number of rejections due to higher trigger latency until the
% analysis will be paused if trig.rejThresh.pause is enabled
trig.rejThresh.note = 5;
% enable to pause the analysis if trig.rejThresh.note is exceeded
trig.rejThresh.pause = 0;
% start trigger of each raw file
trig.startTrig = '254';
% end trigger of each raw file
trig.endTrig = '255';
% response trigger
trig.respTrig = '17';
% miss trigger
trig.missTrig = '99';
% trigger that should be deleted directly after raw data was loaded (due to
% technical issues
trig.delete = {'boundary'};
% trigger codes of events which have to be excluded by definition
% (hypothesis driven) / are already marked as excluded by the stimulation
% scripts 
trig.exclBefore = {'21' '22' '23' '24' '25' '26'};

%% RETRIGGERING

% Change retriggering Matrix in /func/retrigConf.m


%% Triggers & Labels
% 1st Column: Trigger
% 2nd Column: Trigger Label
% 3rd Column: Vector of Tasktypes in which the event occurs
% 4th Column: Color index (same color per condition [across condition comparison])
% 5th Column: Color index (same color per type [within condition comparison])
% 6th Column: Index types belonging together (for choice of same numbers of trials)
% 7th Column: Linestyle ({} = default = '-')
switch lower(method)
  case 'triggers'
    switch taskType
      case 1
        trig.triggers = {        
        % EXAMPLE:
        % Condition 1
        % '120',         {'1stCondTrig1'},           [1 2], [1], [1], [1],{};  
        % Condition 2
        % '210',         {'2ndCondTrig1'},           [1 2], [2], [1], [1],{};
        % '220',         {'2ndCondTrig2'},           [1 2], [2], [2], [1],{};
        % Condition 3
        % '310',         {'3rdCondTrig1'},           [1 2], [3], [1], [1],{};
        % '320',         {'3rdCondTrig2'},           [1 2], [3], [2], [1],{};
        % Condition 4
        % '410',         {'4thCondTrig1'},           [1 2], [4], [1], [1],{};
        % '420',         {'4thCondTrig2'},           [1 2], [4], [2], [1],{};
        % Condition 5
        % {'510' '520'}, {'5thCondBothTrig'},        [1 2], [5], [1], [1],{};
                        };
      case 2
        % ACTIVE TASK (2)
        trig.triggers = {        
        % ... same as above for another task/set of data (preferably
        % active data with response triggers in the EEG)
                        };
    end

    % Define difference waves to compute
    % 1st Column minus 2nd column
    % 3rd Column is name of difference wave
    % 4th Column is color index per condition
    % 5th Column is color index per type
    % 6th Column is line style (default '-')
    trig.diffWaves = {
    % EXAMPLE:
    % '2ndCondTrig1', '2ndCondTrig2', '2ndCondTrig1-minus-2ndCondTrig2', [6], [2], {};
    % '3ndCondTrig1', '3ndCondTrig2', '3ndCondTrig1-minus-3ndCondTrig2', [6], [3], {};
                     };
  case 'retrig'
    %% retrigger events (first of pair/ second of pair)
    trig = {
    %
    %     {'Original trigger'}, {'new trigger'}};  
    %
    % example:
    %     {'41'}, num2str((100 * condOrder(iFile)) + 41);  
             };    
end