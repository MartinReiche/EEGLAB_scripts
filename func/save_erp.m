% combine subject data, calculate amount of usable data and save erpAll
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

function save_erp(subErp,subErpEqual,subTrialInd,corrTrials,numEvent,rejEpochsIn,trialNum,subjects,paths,trig,dur,analysis,filtPar,rejLog,job)

% transform cell array of usable data and erp data to matrices
corrTrials = cell2mat(corrTrials);
subErp = cell2mat(subErp);
subErpEqual = cell2mat(subErpEqual);
% initialize cell array for indices of rejected epochs
rejEpochs = cell(max(subjects),size(trig.triggers,1));

% delete paths.allFiles field
if isfield(paths,'allFiles')
   disp(':: Removing field ''allFiles'' from path structure');
   rmfield(paths,'allFiles');
end

% load chanlocs to save in erp structure
load([paths.resDirAll 'chanlocs.mat']);

% calculate analysis duration
startVec = datevec(dur.Start,'ddd mmm DD HH:MM:SS');
endVec = datevec(dur.End,'ddd mmm DD HH:MM:SS');
dur.secs = etime(endVec,startVec);
dur.s = mod(dur.secs,60);
dur.m = floor(mod((dur.secs/60),60));
dur.h = floor(mod((dur.secs/3600),60));
dur.d = floor(mod((dur.secs/86400),24));
dur.duration = [num2str(dur.d) ' days ' num2str(dur.h) 'h ' num2str(dur.m)...
                'm ' num2str(dur.s) 's'];

% calculate analysis duration for each task if performed on cluster
for iSubj = 1:numel(subjects)
    dur.task(iSubj).subject = subjects(iSubj);
    if analysis.parallel
        startVec = datevec(job.Tasks(iSubj).StartTime,'ddd mmm DD HH:MM:SS');
        endVec = datevec(job.Tasks(iSubj).FinishTime,'ddd mmm DD HH:MM:SS');
        dur.task(iSubj).Start = job.Tasks(iSubj).StartTime;
        dur.task(iSubj).End = job.Tasks(iSubj).FinishTime;
    else
        analysis.parallel
        startVec = datevec(dur.task(iSubj).Start,'ddd mmm DD HH:MM:SS');
        endVec = datevec(dur.task(iSubj).End,'ddd mmm DD HH:MM:SS');
    end        
    dur.task(iSubj).secs = etime(endVec,startVec);
    dur.task(iSubj).s = mod(dur.task(iSubj).secs,60);
    dur.task(iSubj).m = floor(dur.task(iSubj).secs / 60);
    dur.task(iSubj).h = floor(dur.task(iSubj).m / 60);
    dur.task(iSubj).d = floor(dur.task(iSubj).h / 24);
    dur.task(iSubj).duration = [num2str(dur.task(iSubj).d) ' days ' ...
                        num2str(dur.task(iSubj).h) 'h ' num2str(dur.task(iSubj).m) ...
                        'm ' num2str(dur.task(iSubj).s) 's'];
end
%% Initialize ERP matrix
if exist([paths.resDirAll paths.erpFileName],'file')
    disp(' ');
    in = input([':: ERP File already exists, do you want to use this one and ' ...
                'overwrite the data for the given Subjects? (Y/n)'],'s');
    answ = 1;
    while answ
        if (strcmpi(in,'Y')) || isempty(in)
            disp(':: Loading ERP file');
            load([paths.resDirAll paths.erpFileName]);
            % retrieve old erpAll and rejepoch matrices (new subject data
            % will be inserted in the respective line)
            erpAll = erp.erpAll;
            erpAllEqual = erp.erpAllEqual;
            rejEpochs = erp.rejEpochs;
            answ = 0;
            % save analysis duration in minutes, subjects, date of analysis, machine and
            % EEGLAB version
            erp.meta(end+1).subjects = subjects;
            erp.meta(end).startTime = dur.Start;
            erp.meta(end).finishTime = dur.End;
            erp.meta(end).duration = dur.duration;
            erp.meta(end).taskDur = dur.task;
            erp.meta(end).version = ['EEGLAB ' eeg_getversion];
            erp.meta(end).machine = system_dependent('getos');
            % store configuration parameters of current analysis (append
            % and store separately for each analysis)
            erp.analysis(end+1) = analysis;
            erp.paths(end+1) = paths;
            erp.trig(end+1) = trig;
            erp.filtPar(end+1) = filtPar;            
        elseif strcmpi(in,'n')
            disp(':: Creating new ERP file');
            % initialize erpAll matrices
            erpAll = zeros(max(subjects),...
                           size(trig.triggers,1),...
                           size(subErp,3),...
                           size(subErp,4));
            erpAllEqual = erpAll;
            % initialize cell array for indices of rejected epochs
            rejEpochs = cell(max(subjects),size(trig.triggers,1));
            % initialize trial indices for equal amount of trials per
            % trigger group (for each subject and each trigger)
            erp.subTrialInd = cell(max(subjects),1);
            % initialize number of correct trials per subject and trigger
            erp.corrTrials = cell(max(subjects),1);
            % initialize number of overall trials per subject and trigger
            erp.eventCount = cell(max(subjects),1);
            % initialize number of trials per trigger group for
            % each subject
            erp.trialNum = cell(max(subjects),1);
            answ = 0;
            % save meta data
            erp.meta.subjects = subjects;
            erp.meta(end).startTime = dur.Start;
            erp.meta(end).finishTime = dur.End;
            erp.meta(end).duration = dur.duration;
            erp.meta(end).taskDur = dur.task;
            erp.meta.version = ['EEGLAB ' eeg_getversion];
            erp.meta.machine = system_dependent('getos');
            % store configuration parameters of current analysis
            erp.analysis = analysis;
            erp.paths = paths;
            erp.trig = trig;
            erp.filtPar = filtPar;           
        else
            in = input([':: ' in ' is not an option! (Y/n)'],'s');
        end
    end
else
    
    disp(':: Did not find file, creating new one.');
    % initialize erpAll matrices
    erpAll = zeros(max(subjects),...
                   size(trig.triggers,1),...
                   size(subErp,3),...
                   size(subErp,4));
    erpAllEqual = erpAll;
    % initialize cell array for indices of rejected epochs
    rejEpochs = cell(max(subjects),size(trig.triggers,1));
    % initialize trial indices for equal amount of trials per
    % trigger group (for each subject and each trigger)
    erp.subTrialInd = cell(max(subjects),1);
    % initialize number of correct trials per subject and trigger
    erp.corrTrials = cell(max(subjects),1);
    % initialize number of overall trials per subject and trigger
    erp.eventCount = cell(max(subjects),1);
    % initialize number of trials per trigger group for
    % each subject
    erp.trialNum = cell(max(subjects),1);
    % save meta data
    erp.meta.subjects = subjects;
    erp.meta(end).startTime = dur.Start;
    erp.meta(end).finishTime = dur.End;
    erp.meta(end).duration = dur.duration;
    erp.meta(end).taskDur = dur.task;
    erp.meta.version = ['EEGLAB ' eeg_getversion];
    erp.meta.machine = system_dependent('getos');
    % store configuration parameters of current analysis
    erp.analysis = analysis;
    erp.paths = paths;
    erp.trig = trig;
    erp.filtPar = filtPar;            
end

%% Combine Subjects and get usable data
% initialize folder for rejection information
paths.rejFolder = [paths.resDirAll 'usable_Data_' paths.erpFileName(1:end-4) '/'];

if ~exist(paths.rejFolder,'dir');
    disp(':: Creating new folder to store usable Data text files');
    mkdir(paths.rejFolder)
end

% open text file for summary of rejection
rejSum = fopen([paths.rejFolder 'Rejected_Data_Log.txt'],'w');
fprintf(rejSum,'%s\n','-------------------------------------------------------------------------------------------');
fprintf(rejSum,'%s\n\n\n',['                          :: REJECTION LOG FILE FOR ALL SUBJECTS ::']);
% add meta information of preprocessing to rejection log file header
fprintf(rejSum,'%s\n',['Start of preprocessing:       ' dur.Start]);
fprintf(rejSum,'%s\n',['End of preprocessing:         ' dur.End]);
fprintf(rejSum,'%s\n\n',['Duration of preprocessing:    ' dur.duration]);
if analysis.parallel
    fprintf(rejSum,'%s\n',['Cluster configuration:        ' analysis.core]);
else
    fprintf(rejSum,'%s\n',['Cluster configuration:        ' 'serial']);
end
fprintf(rejSum,'%s\n',['EEGLAB version:               ' erp.meta(end).version]);
fprintf(rejSum,'%s\n',['Machine:                      ' erp.meta(end).machine]);
fprintf(rejSum,'%s\n','-------------------------------------------------------------------------------------------');
fprintf(rejSum,'%s\n',' ');

disp(':: Combining subject data');
% run through all Subjects
for iSubj = 1:size(subjects,2)

    % get number of events for current Subjects
    eventCount = numEvent{iSubj};
    
    %% Assign erpSub data to erpAll and fetch indices of rejected epochs %%%%%%%%%%%%%%%%%%%%%%
    erpAll(subjects(iSubj),:,:,:) = subErp(iSubj,:,:,:);
    erpAllEqual(subjects(iSubj),:,:,:) = subErpEqual(iSubj,:,:,:);
    rejEpochs(subjects(iSubj),:) = rejEpochsIn(iSubj,:);
    erp.corrTrials{subjects(iSubj),1} = corrTrials(iSubj,:);
    erp.subTrialInd{subjects(iSubj),1} = subTrialInd{iSubj,1};
    erp.eventCount{subjects(iSubj),1} = numEvent{iSubj,1};
    erp.trialNum{subjects(iSubj),1} = trialNum{iSubj,1};
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Report usable data per trigger
    % initialize folder for rejection information
    paths.rejFolder = [paths.resDirAll 'usable_Data_' paths.erpFileName(1:end-4) '/'];
    
    if ~exist(paths.rejFolder,'dir');
        disp(':: Creating new folder to store usable Data text files');
        mkdir(paths.rejFolder);
    end
    % open text file for report on rejected epoches
    rejFile = fopen([paths.rejFolder 'Subj' num2str(subjects(iSubj), '%0.2d')...
                     'reject.txt'],'w');

    % write rejection file heading (containing labels of predefined excluded channels)

    fprintf(rejFile,'%s\n','-------------------------------------------------------------------------------------------');
    fprintf(rejFile,'%s\n\n\n',['                          :: REJECTION LOG FILE OF SUBJECT ' num2str(subjects(iSubj), '%0.2d') ' ::']);
    
    if ~isempty(rejLog(iSubj).exclVector)
        exclLabels = [];
        % get labels of excluded electrodes
        for iChannel = rejLog(iSubj).exclVector
            exclLabels = [exclLabels ' ' chanlocs(iChannel).labels];
        end
        fprintf(rejFile,'%s\n',['Rejection mode:               ' analysis.rejLabel{analysis.rejmode+1,:} ', excluding channels' exclLabels]);
    else
        fprintf(rejFile,'%s\n\n',['Rejection mode:               ' analysis.rejLabel{analysis.rejmode+1,:} ', including all channels']);
    end

    % add meta information of preprocessing to rejection log file header
    fprintf(rejFile,'%s\n',['Start of preprocessing:       ' dur.Start]);
    fprintf(rejFile,'%s\n',['End of preprocessing:         ' dur.End]);
    fprintf(rejFile,'%s\n\n',['Duration of preprocessing:    ' dur.duration]);
    fprintf(rejFile,'%s\n',['Start of task:                ' dur.task(iSubj).Start]);
    fprintf(rejFile,'%s\n',['End of task:                  ' dur.task(iSubj).End]);
    fprintf(rejFile,'%s\n\n',['Duration of task:             ' dur.task(iSubj).duration]);
    % add cluster and machine configuration
    if analysis.parallel
        fprintf(rejFile,'%s\n',['Cluster configuration:        ' analysis.core]);
    else
        fprintf(rejFile,'%s\n',['Cluster configuration:        ' 'serial']);
    end
    fprintf(rejFile,'%s\n',['EEGLAB version:               ' erp.meta(end).version]);
    fprintf(rejFile,'%s\n\n',['Machine:                      ' erp.meta(end).machine]);
    
    if analysis.rejFlatepochs
        fprintf(rejFile,'%s\n',['Channels with rejection proportions exceeding ' num2str(analysis.chanMaxRej*100) ' % of epoches are reported separately']);
    end
    fprintf(rejFile,'%s\n\n','-------------------------------------------------------------------------------------------');

    for iTrig = 1:size(trig.triggers,1)
        
        % convert cell of trigger names to concatenated trigger string in
        % case several triggers are selected at once
        if iscell(trig.triggers{iTrig,1})
            currTrig = num2str(cell2mat(trig.triggers{iTrig,1}),'%0.2d');
        else
            currTrig = num2str(trig.triggers{iTrig,1},'%0.2d');
        end
                
        % Get Proportion of usable Data for current subject and trigger
        corrProp(iSubj,iTrig) = corrTrials(iSubj,iTrig)/eventCount(iTrig);
        if isnan(corrProp(iSubj,iTrig))
            corrProp(iSubj,iTrig) = 0; 
        end
        
    end

    % report overall amount of usable data
    fprintf(rejSum,'%s\n',[':: Subject ' num2str(subjects(iSubj), '%0.2d') ': ' num2str(round(mean(corrProp(iSubj,:))*1000)/10) ' %']);
    fprintf(rejFile,'%s\n\n',['USABLE DATA OF SUBJECT ' num2str(subjects(iSubj), '%0.2d') ': ' num2str(round(mean(corrProp(iSubj,:))*1000)/10)...
                        '%']);
    fprintf(rejFile,'%s\n\n','-------------------------------------------------------------------------------------------');


    
    for iTrig = 1:size(trig.triggers,1)

        % convert cell of trigger names to concatenated trigger string in
        % case several triggers are selected at once
        if iscell(trig.triggers{iTrig,1})
            currTrig = num2str(cell2mat(trig.triggers{iTrig,1}),'%0.2d');
        else
            currTrig = num2str(trig.triggers{iTrig,1},'%0.2d');
        end

        
        % Display usable Data for current subject and trigger and store in rejFile
        disp(['Trigger ' currTrig ' ' num2str(corrTrials(iSubj,iTrig)) ' Events averaged = '...
              num2str(round(corrProp(iSubj,iTrig)*100)) '%']);
        fprintf(rejFile,'%s\n',['Trigger ' currTrig ' Events averaged:       '...
                            num2str(round(corrProp(iSubj,iTrig)*100)) '%']);
        fprintf(rejFile,'%s\n',['  - number of epoches (original):   ' num2str(eventCount(iTrig))]);
        fprintf(rejFile,'%s\n',['  - number of epoches   (usable):   ' num2str(corrTrials(iSubj,iTrig))]);
        
        % look up whether epoches were rejected for being flat for the
        % current trigger of the current subject
        for iLine = 1:size(rejLog(iSubj).flat,1)
            if rejLog(iSubj).flat{iLine,1} == iSubj && rejLog(iSubj).flat{iLine,2} == iTrig && rejLog(iSubj).flat{iLine,3}
                % and if so, report the number of epoches in the rejection
                % log file
                
                % report amount of averaged trials
                fprintf(rejFile,'%s\n',['  - theoretically ' num2str(rejLog(iSubj).flat{iLine,4}) ' trials, practically ' ...
                                    num2str(rejLog(iSubj).flat{iLine,5}) ' trials']);
                if analysis.rejFlatepoches
                    fprintf(rejFile,'%s\n',['  -- ' num2str(rejLog(iSubj).flat{iLine,3}) ' trials rejected for being flat']);
                end
                
            end
        end

        % report minimum and maximum amplitude changes for trigger
        fprintf(rejFile,'%s\n',['  -- minimum amplitude change:      ' num2str(rejLog(iSubj).minChange(iTrig)) ' microV']);
        fprintf(rejFile,'%s\n',['  -- maximum amplitude change:      ' num2str(rejLog(iSubj).maxChange(iTrig)) ' microV']);
        
        % look for specific channels causing rejection for current trigger
        for iLine = 1:size(rejLog(iSubj).chanReport,1)
            if strcmpi(rejLog(iSubj).chanReport(iLine,2),currTrig)
                % if channels were detected for causing rejection of the given trigger, report
                % them in the rejection report file
                fprintf(rejFile,'%s\n',['  --- channel ' num2str(rejLog(iSubj).chanReport{iLine,3}) ' (' rejLog(iSubj).chanReport{iLine,4} ') '...
                                    'rejected Epochs: ' num2str(round(rejLog(iSubj).chanReport{iLine,5}*100)) ...
                                    ' %']);
            end
        end
        
        % add empty line after specific channel rejection log
        fprintf(rejFile,'%s\n','');

        disp(' ');
    end
    % close rejection log file of current subject
    fclose(rejFile);
    
    % save vectors of rejected epochs for all subjects
    disp(':: Saving rejected epoch indices.');

end

% display and store average of usable data across subjects and conditions
for iSubj = 1:numel(subjects)
    disp(' ');
    disp([':: Subject ' num2str(subjects(iSubj), '%0.2d') ' ' num2str(round(mean(corrProp(iSubj,:))*1000)/10) '% usable data across conditions ']);
    for iTrig = 1:size(trig.triggers,1)
        
        % convert cell of trigger names to concatenated trigger string in
        % case several triggers are selected at once
        if iscell(trig.triggers{iTrig,1})
            currTrig = num2str(cell2mat(trig.triggers{iTrig,1}),'%0.2d');
        else
            currTrig = num2str(trig.triggers{iTrig,1},'%0.2d');
        end
        
        if corrProp(iSubj,iTrig) < 0.8
            disp(['subject ' num2str(iSubj, '%0.2d') ' trigger ' currTrig ' '...
                  num2str(corrTrials(iSubj,iTrig)) ' trials averaged = ' num2str(round(corrProp(iSubj,iTrig)*100)) '%']);
        end
    end
end
% write overall amount of usable data to rejection summary file
fprintf(rejSum,'%s\n',' ');
fprintf(rejSum,'%s\n','-------------------------------------------------------------------------------------------');
fprintf(rejSum,'%s\n',['USABLE DATA ACROSS SUBJECTS AND CONDITIONS:       ' ...
                    num2str(round(mean(mean(corrProp,2))*1000)/10) ' % ' ...
                    '(based on condition means)']);
fprintf(rejSum,'%s\n',['USABLE DATA ACROSS SUBJECTS AND CONDITIONS:       ' ...
                    num2str(round(sum(sum(corrTrials))/(size(subjects,2)* ...
                                                  sum(eventCount))*1000)/10) ...
                    ' % (based on trial numbers overall)']);
fprintf(rejSum,'%s\n','-------------------------------------------------------------------------------------------');

% close rejection summary file
fclose(rejSum);
% display stats of usable data over all subjects
disp(' ');
disp('based on condition means');
disp([num2str(round(mean(mean(corrProp,2))*1000)/10) '% usable data across subjects and conditions ']);
disp('based on trial numbers overall');
disp([num2str(round(sum(sum(corrTrials))/(size(subjects,2)*sum(eventCount))*1000)/10) '% usable data across subjects and conditions ']);

% save erpAll
erp.erpAll = erpAll;
erp.erpAllEqual = erpAllEqual;
% save rejected epoch indices
erp.rejEpochs = rejEpochs;
% save chanlocs in erp structure
erp.chanlocs = chanlocs;
% save erp structure
disp(':: Saving results, this might take some time.');
save([paths.resDirAll paths.erpFileName],'erp');
% Prompt to wait for User
input('continue with enter');
