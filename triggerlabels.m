% Names of triggers and trigger configurations
function trig = triggerlabels(trig,taskType)
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
trig.startTrig = 'S254';
% end trigger of each raw file
trig.endTrig = 'S255';
% response trigger
trig.respTrig = ' ';
% miss trigger
trig.missTrig = '99';
% trigger codes of events which have to be excluded by definition
% (hypothesis driven) / are already marked as excluded by the stimulation
% scripts 
trig.exclBefore = [ ];

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

switch taskType
  case 1
    % PASSIVE TASK (1)
    trig.triggers = {
    % condition 1 (1st & 2nd tone of "pair")
        'S 11',        {'first-tone-1'},          [1]  , [1], [1], [1], {};
        'S 12',        {'second-tone-1'},         [1]  , [2], [1], [1], {};
    % condition 1 (omission of 1st & 2nd tone of "pair")
        'S 13',        {'first-omission-1'},      [1]  , [3], [1], [2], {};
        'S 14',        {'second-omission-1'},     [1]  , [4], [1], [2], {};
    % condition 2 (1st & 2nd tone of "pair")
        'S 21',        {'first-tone-2'},          [1]  , [1], [2], [1], {};
        'S 22',        {'second-tone-2'},         [1]  , [2], [2], [1], {};
    % condition 2 (omission of 1st & 2nd tone of "pair")
        'S 23',        {'first-omission-2'},      [1]  , [3], [2], [2], {};
        'S 24',        {'second-omission-2'},     [1]  , [4], [2], [2], {};
    % condition 3 (1st & 2nd tone of "pair")
        'S 31',        {'first-tone-3'},          [1]  , [1], [3], [1], {};
        'S 32',        {'second-tone-3'},         [1]  , [2], [3], [1], {};
    % condition 3 (omission of 1st & 2nd tone of "pair")
        'S 33',        {'first-omission-3'},      [1]  , [3], [3], [2], {};
        'S 34',        {'second-omission-3'},     [1]  , [4], [3], [2], {};
    % condition 4 (1st & 2nd tone of "pair")
        'S 41',        {'first-tone-4'},          [1]  , [1], [4], [1], {};
        'S 42',        {'second-tone-4'},         [1]  , [2], [4], [1], {};
    % condition 4 (omission of 1st & 2nd tone of "pair")
        'S 43',        {'first-omission-4'},      [1]  , [3], [4], [2], {};
        'S 44',        {'second-omission-4'},     [1]  , [4], [4], [2], {};
    % condition 5 (1st & 2nd tone of "pair")
        'S 51',        {'first-tone-5'},          [1]  , [1], [5], [1], {};
        'S 52',        {'second-tone-5'},         [1]  , [2], [5], [1], {};
    % condition 5 (omission of 1st & 2nd tone of "pair")
        'S 53',        {'first-omission-5'},      [1]  , [3], [5], [2], {};
        'S 54',        {'second-omission-5'},     [1]  , [4], [5], [2], {};
                    };
  case 2
        trig.triggers = {
    % ACTIVE TASK (2)
    % THERE IS NO ACTIVE TASK FOR THIS EXPERIMENT
                        };
        end

% Define difference waves to compute
% 1st Column minus 2nd column
% 3rd Column is name of difference wave
% 4th Column is color index per condition
% 5th Column is color index per type
% 6th Column is line style (default '-')

trig.diffWaves = {
    'first-tone-1', 'second-tone-1', 'tone-diff-1', [5], [1], {};
    'first-tone-2', 'second-tone-2', 'tone-diff-2', [5], [2], {};
    'first-tone-3', 'second-tone-3', 'tone-diff-3', [5], [3], {};
    'first-tone-4', 'second-tone-4', 'tone-diff-4', [5], [4], {};
    'first-tone-5', 'second-tone-5', 'tone-diff-5', [5], [5], {};
                 };

