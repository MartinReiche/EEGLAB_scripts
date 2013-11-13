% Interpolate single channels depending on block and subject 

% Author: Martin Reiche, martin.reiche@uni-oldnburg.de
% Based on a script by Alexandra Bendixen

function EEG = interpChan(EEG,analysis,nSubj,iFile)
% go through the chaninterp array
if ~isempty(analysis.chanInterp)
    for iInterp = 1:size(analysis.chanInterp,1)
        % if the current subject and block is in the current line of the analysis.chanInterp
        % array and the

        if ~ischar(analysis.chanInterp{iInterp,2})
            analysis.chanInterp{iInterp,2} = EEG.chanlocs(analysis.chanInterp{iInterp,2}).labels;
        end

        if (nSubj == analysis.chanInterp{iInterp,1}) && any(ismember(iFile, ...
                                                         analysis.chanInterp{iInterp,3}))

            % find the channel number for the channel to be interpolated
            for iChan = 1:size(EEG.chanlocs,2)

                if strcmpi(EEG.chanlocs(iChan).labels,analysis.chanInterp{iInterp,2})
                    disp(' ');
                    disp([':: Interpolating channel ' analysis.chanInterp{iInterp,2} ' ' ...
                          ' for Subject '  num2str(nSubj,'%0.2d') ...
                          ' Block ' num2str(iFile,'%0.2d')]);
                    disp(' ');
                    interpolLabel = EEG.chanlocs(iChan).labels;
                    interpolNum = iChan;
                    
                    hpos = []; vpos = [];
                    %find eye channels
                    for n=1:size(EEG.chanlocs,2)
                        if strcmp(EEG.chanlocs(n).labels,'HEOG')
                            hpos = n;
                        end
                        if strcmp(EEG.chanlocs(n).labels,'VEOG')
                            vpos = n;
                        end
                    end
                    if isempty(hpos)
                        error('failed to find HEOG');
                    end
                    if isempty(vpos)
                        error('failed to find VEOG');
                    end
                    EEG2 = pop_select(EEG,'nochannel',[hpos vpos]);

                    %look for position in new EEG2 structure
                    eFound = 0;
                    for iChannel = 1:size(EEG2.chanlocs,2) %determine relevant channel numbers
                        if strcmp(EEG2.chanlocs(iChannel).labels,interpolLabel)
                            eFound = 1;
                            interpolNum2 = iChannel;
                        end
                    end
                    if ~eFound
                        error(['failed to find ' char(interpolLabel)]);
                    end

                    % interpolate
                    EEG2 = eeg_interp(EEG2,interpolNum2);
                    disp(' ');
                    disp([':: Replacing data of electrode ' ...
                          EEG.chanlocs(interpolNum).labels ' with interpolated data from ' EEG2.chanlocs(interpolNum2).labels]);
                    disp(' ');
                    EEG.data(interpolNum,:,:) = EEG2.data(interpolNum2,:,:);

                    
                end
            end
        end
    end


end

end