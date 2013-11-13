% Merge the Data Files of given Subject
%
% Copyright (c) 2013 Martin Reiche, Carl-von-Ossietzky-University Oldenburg
% Author: Martin Reiche, martin.reiche@uni-oldnburg.de

function merge_data(iSubj,paths,trig,eventCount,condOrder)
    
%% Merge all files belonging to one block in one file
    
% delete all epoch files which contain only one trial
    allFiles = dir([paths.resDir paths.resFileSubSpec '*'...
                    paths.resFileBlockSpec '*'...
                    paths.resFileTrigSpec '*'...
                    paths.resFileExt]);
    for iFile = 1:size(allFiles,1)
        EEG = pop_loadset(allFiles(iFile).name,paths.resDir);
        if EEG.trials == 1
            disp(' ');
            disp(' ');
            disp([':: Only one epoch detected in ' allFiles(iFile).name]);
            disp([':: Cannot merge this file, deleting it']);
            disp(' ');
            
            % if there is only one epoch in a filpe, automatically
            % reject this epoch and save the filename 
            
            fid = fopen([paths.resDir paths.resFileSubSpec num2str(iSubj, '%0.2d') ...
                         'single_epoch.txt'],'a');
            fprintf(fid, '%s\n', allFiles(iFile).name);
            % close single epoch rejection file
            fclose(fid);
            delete([paths.resDir allFiles(iFile).name(1:end-length(paths.resFileExt)) '*']);
        end
    end

    % start merging
    disp(':: Merging Trigger Files');
    for iTrig = 1:size(trig.triggers,1)
        % for all triggers
        trigFiles = dir([paths.resDir paths.resFileSubSpec '*'...
                         paths.resFileBlockSpec '*'...
                         paths.resFileTrigSpec...
                         num2str(trig.triggers{iTrig,1},'%0.2d')...
                         paths.resFileExt]);
        
        
        % create list of EEG sets
        setList = 1:size(trigFiles,1);
        
        for iFile = 1:size(trigFiles,1)
            ALLEEG(iFile) = pop_loadset(trigFiles(iFile).name,paths.resDir);
        end
    
        % if there were epoch files for current trigger in current condition
        if exist('ALLEEG','var') 
            if size(ALLEEG,2) > 1
                EEG = pop_mergeset(ALLEEG,setList,0);
            else
                EEG = ALLEEG(1);
            end
            EEG = pop_saveset(EEG,[paths.resFileSubSpec  num2str(iSubj,'%0.2d') ...
                                paths.resFileTrigSpec...
                                num2str(trig.triggers{iTrig,1},'%0.2d') ...
                                'all.set'],paths.resDir);
            clear ALLCOM ALLEEG CURRENTSET EEG LASTCOM;
        end
        
    end
    
    % remove old files
    disp(' ');
    disp(':: Delete old files');
    delete([paths.resDir paths.resFileSubSpec num2str(iSubj,'%0.2d')... 
            paths.resFileBlockSpec '*']);
    
    % compute initial number of all events per trigger
    eventCount = sum(eventCount,2);
    
    for iCond = 1:max(condOrder)
        numEvent(iCond,:) = sum(eventCount(eventCount(:,end) == iCond,1:end-1),1);
    end
    % save numbers of events in result dir of current subject
    disp(':: Saving initial number of events');
    save([paths.resDir paths.resFileSubSpec num2str(iSubj, '%0.2d') 'nTrials.mat'],'eventCount');
    
end