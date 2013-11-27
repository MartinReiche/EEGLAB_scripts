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

function save_erp(subErp,subErpEqual,subTrialInd,corrTrials,numEvent,rejEpochsIn,trialNum,subjects,paths,trig,tStart,tEnd,analysis,filtPar)

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

%% Initialize ERP matrix
if ~analysis.batchMode && exist([paths.resDirAll paths.erpFileName],'file')
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
            erp.meta(end).duration = (tEnd/60);
            erp.meta(end).date = datestr(now);
            erp.meta(end).version = ['EEGLAB ' eeg_getversion];
            erp.meta(end).maschine = system_dependent('getos');
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
            erp.meta.duration = (tEnd/60);
            erp.meta.date = datestr(now);
            erp.meta.version = ['EEGLAB ' eeg_getversion];
            erp.meta.maschine = system_dependent('getos');
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
    
    paths

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
    erp.meta.duration = (tEnd/60);
    erp.meta.date = datestr(now);
    erp.meta.version = ['EEGLAB ' eeg_getversion];
    erp.meta.maschine = system_dependent('getos');
    % store configuration parameters of current analysis
    erp.analysis = analysis;
    erp.paths = paths;
    erp.trig = trig;
    erp.filtPar = filtPar;            
end

%% Combine Subjects and get usable data
% initialize folder for rejection information
if analysis.batchMode
    % rename ERP File
    paths.rejFolder = [paths.resDirAll 'JOB_' num2str(analysis.jobIndex) '_usable_Data_' paths.erpFileName(1:end-4) '/'];
    paths.erpFileName = ['JOB_' num2str(analysis.jobIndex) '_'  paths.erpFileName]; 
else
    paths.rejFolder = [paths.resDirAll 'usable_Data_' paths.erpFileName(1:end-4) '/'];
end
if ~exist(paths.rejFolder,'dir');
    disp(':: Creating new folder to store usable Data text files');
    mkdir(paths.rejFolder)
end

disp(':: Combining subject data');
% run through all Subjects
for iSubj = 1:size(subjects,2)

    % open text file for report on rejected epoches
    rejFile = fopen([paths.rejFolder 'Subj' num2str(subjects(iSubj), '%0.2d')...
                     'reject.txt'],'a');
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
        
        % Display usable Data for current subject and trigger and store in rejFile
        disp(['Trigger ' currTrig ' ' num2str(corrTrials(iSubj,iTrig)) ' Events averaged = '...
              num2str(round(corrProp(iSubj,iTrig)*100)) '%']);
        fprintf(rejFile,'%s\n',['Trigger ' currTrig ' ' currTrig ' Events averaged = '...
                            num2str(round(corrProp(iSubj,iTrig)*100)) '%']);
        disp(' ');
    end
    
    fprintf(rejFile,'%s\n','----------------------------------------------------------------');
    fprintf(rejFile,'%s\n\n',[':: Subject ' num2str(subjects(iSubj), '%0.2d') ' ' num2str(round(mean(corrProp(iSubj,:))*1000)/10)...
                        '% usable data across conditions ']);  
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
disp(' ');
disp('based on condition means');
disp([num2str(round(mean(mean(corrProp,2))*1000)/10) '% usable data across subjects and conditions ']);
disp('based on trial numbers overall');
disp([num2str(round(sum(sum(corrTrials))/(size(subjects,2)*sum(eventCount))*1000)/10) '% usable data across subjects and conditions ']);

% load chanlocs to save in erp structure
load([paths.resDirAll 'chanlocs.mat']);
% save erpAll
erp.erpAll = erpAll;
erp.erpAllEqual = erpAllEqual;
% save rejected epoch indices
erp.rejEpochs = rejEpochs;
% save chanlocs in erp structure
erp.chanlocs = chanlocs;
% save erp structure
disp(':: Saving results');
save([paths.resDirAll paths.erpFileName],'erp');
% remove temporarily saved chanlocs
delete([paths.resDirAll 'chanlocs.mat']);
% Prompt to wait for User
if ~analysis.batchMode
    input('continue with enter');
end