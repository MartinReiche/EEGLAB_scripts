% All the plotting of ERP curves (specified in config.m) is handled in
% this function. Plotting is handled in two stages: first plots of
% specified curves on all (or a specified set of) electrodes will be
% presented per figure. In the socond step, curves will be presented on a
% specified electrode and stats will be shown for specified windows.
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

function plot_erp(erpAll,chanlocs,plotPar,trig,analysis,paths,taskType,restoredConf)

% use saved erp window parameters
    if analysis.savedErpWin && (nargin > 7)
        plotPar.xScale = restoredConf.analysis(end).erpWin;
        analysis.erpWin = restoredConf.analysis(end).erpWin;
        analysis.baseWin = analysis.erpWin(1);
    end
    
    % load current trigger parameters
    currTrigPars = config('triggers','task', taskType);

    % Save input arguments before manipulating
    urerpAll = erpAll;
    urchanlocs = chanlocs;
    urtrig = trig;
    
    % 1st first take the labels without diference wave and form a second column
    % and numerate from 1 to end (should correspond with the indices in
    % erpAll). Take trig.diffInd which should hold the erpAll Indices and
    % vertically concetenate it with the one above.

    % if there are no difference waves
    if isempty(trig.diffInd)
        trig.trigLabels = {trig.trigLabels, (1:size(trig.trigLabels,1))'};
        labels = trig.trigLabels{1};
        indices = trig.trigLabels{2};
        trig.color = trig.color;
    else
        trig.trigLabels = {trig.trigLabels, (1:size(trig.trigLabels,1))'};
        labels = [trig.trigLabels{1};trig.diffInd(:,2)];
        indices = [trig.trigLabels{2}; cell2mat(trig.diffInd(:,3))];
        trig.color = [trig.color; trig.diffCol];
    end

    %% get the line style of the curves from current config file
    plotPar.style = {};
    for iCurve = 1:size(labels,1)
        foundLabel = 0;
        % search trigger array for label
        for iTrig = 1:size(currTrigPars.triggers,1)
            if strcmpi(currTrigPars.triggers{iTrig,2},labels(iCurve))
                foundLabel = 1;
                if isempty(currTrigPars.triggers{iTrig,7})
                    plotPar.style{iCurve,1} = '-';
                else
                    plotPar.style(iCurve,1) = currTrigPars.triggers{iTrig,7};
                end
            end
        end
        % search diffwaves array for label
        for iTrig = 1:size(currTrigPars.diffWaves,1)
            if strcmpi(currTrigPars.diffWaves{iTrig,3},labels(iCurve))
                foundLabel = 1;
                if isempty(currTrigPars.diffWaves{iTrig,6})
                    plotPar.style{iCurve,1} = '-';
                else
                    plotPar.style(iCurve,1) = currTrigPars.diffWaves{iTrig,6};
                end
            end
        end
        if ~foundLabel
            plotPar.style{iCurve,1} = '-';
        end
    end

    % compute grand average ERPs
    disp(':: Compute grand average ERPs');
    gavr = mean(erpAll,1);

    % check validity of component windows
    if ~isempty(plotPar.compWin)
        for nComp = size(plotPar.compWin,1):-1:1
            if analysis.erpWin(1) > plotPar.compWin(nComp,1) || analysis.erpWin(2) < plotPar.compWin(nComp,2)
                disp(':: WARNING: Component window exceeds epoch range, skipping.');
                plotPar.comps(nComp) = [];
                plotPar.compWin(nComp,:) = [];
            end
        end
    else
        disp(':: WARNING: No component window given, skipping statistics');
    end
    if isempty(plotPar.compWin)
        plotPar.comps = {};
        plotPar.compWin = [];
    end
    
    if analysis.plotERPflag
        %% Plot curves at specified set of channels per figure
        % get position of channels for plotting
        channels2plot = zeros(1,size(plotPar.plotChannels,1));
        for nChannel = 1:size(plotPar.plotChannels,1) 
            for n=1:numel(chanlocs)
                if strcmp(chanlocs(n).labels,plotPar.plotChannels(nChannel,:))
                    channels2plot(nChannel) = n;
                end
            end
        end
        % determine which cells in plotPar.plotConds belong to the current task
        % and rearrange it
        plotConds = {};
        for iCond = 1:size(plotPar.plotConds,2)
            % go through each cell in plotPar.plotConds
            if any(ismember(taskType,cell2mat(plotPar.plotConds{iCond}(end))))
                plotConds{end+1} = plotPar.plotConds{iCond}(1:end-1); 
            end
        end

        % get colors
        if plotPar.autoColor
            % plotColor = distinguishable_colors(size(plotConds{iCond},2));
            plotColor = distinct_color(size(plotConds{iCond},plotPar.plotCondsCol));
        else
            % plotColor = distinguishable_colors(max(trig.color(:,2)));
            plotColor = distinct_color(max(trig.color(:,plotPar.plotCondsCol)));
        end
        % get the subplot parameters (dimensions and iterators)
        plotDim = getDisp('electrode array','parameters',plotPar,'structure',channels2plot,'legend',1);    

        % plot one figure for each cell in plotPar.plotConds (each condition) 
        for iCond = 1:size(plotConds,2)
            % go through each cell of plotConds
            fh = figure('WindowButtonDownFcn',@closeFig);
            plotPar.figInd.curr = iCond;
            % initialize index variable
            currInd = [];
            plotPar.currColor = [];
            plotPar.currStyle = [];
            % get the data to plot 
            for iCurve = 1:size(plotConds{iCond},2) 
                % go through each Curve of the current figure
                foundLabel = 0;
                for iIndex = 1:size(labels,1)
                    % check current curve label against each label of the
                    % curves stored in erpAll to get the index of the current
                    % curve label 
                    if strcmp(plotConds{iCond}(iCurve),labels(iIndex))
                        foundLabel = 1;
                        currInd = [currInd iIndex];
                        if plotPar.autoColor
                            plotPar.currColor = [plotPar.currColor; plotColor(iCurve,:)];
                        else
                            plotPar.currColor = [plotPar.currColor; plotColor(trig.color(iIndex,plotPar.plotCondsCol),:)];
                            plotPar.currStyle = [plotPar.currStyle; plotPar.style(iIndex)];
                        end
                    end        
                end
                if isempty(currInd)
                    % through error if requested curve label is not availabel
                    error([':: Did not find curve label ' plotConds{iCond}{iCurve}]);
                end
            end
            % get the labels of the current plot
            currLabels = plotConds{iCond};
            
            % get minimal and maximal value for current figure
            iChan = [];
            maxVal = getMax('electrode array',erpAll,analysis,plotPar,plotConds{iCond},labels,chanlocs,iChan,channels2plot);
            plotPar.yScale(1) = maxVal.erpMin;
            plotPar.yScale(2) = maxVal.erpMax;
            plotPar.yScaleMean(1) = maxVal.meanMin;
            plotPar.yScaleMean(2) = maxVal.meanMax;

            for iChan = 1:size(channels2plot,2)
                plotPar.currChan = iChan;
                plotPar.currChanLabel = plotPar.plotChannels{iChan};
                chanData = [];            
                spData.resAll = [];
                spData.resHead = {};
                spData.erpMean = zeros(size(plotPar.compWin,1),size(currInd,2));
                spData.erpErr = spData.erpMean;
                % get the data for the current condition and the current channel
                for iCurve = 1:size(currInd,2)
                    chanData(:,iCurve) = squeeze(gavr(:,currInd(iCurve),:,channels2plot(iChan)));
                    % Get Data for saving
                    for iWin = 1:size(plotPar.compWin,1)
                        % get the current window
                        statWin = [round(((plotPar.compWin(iWin,1)+abs(plotPar.xScale(1)))*analysis.sampRate)/1000) ...
                                   round(((plotPar.compWin(iWin,2)+abs(plotPar.xScale(1)))*analysis.sampRate)/1000)];
                        % Save Data in results matrix columns resemble conditions 
                        spData.resAll = [spData.resAll mean(squeeze(erpAll(:,currInd(iCurve),statWin(1):statWin(2),channels2plot(iChan))),2)];
                        % write header for current column
                        spData.resHead = [spData.resHead [char(plotPar.plotConds{iCond}(iCurve)) 'Win' num2str(iWin)]];
                        spData.erpMean(iWin,iCurve) = mean(mean(squeeze(erpAll(:,currInd(iCurve),statWin(1):statWin(2),channels2plot(iChan))),2));
                        spData.erpErr(iWin,iCurve) = std(mean(squeeze(erpAll(:,currInd(iCurve),statWin(1):statWin(2),channels2plot(iChan))),2))/sqrt(size(erpAll,1));
                    end 
                end
                % create the subplot for the current channel
                erpSpHandle = subplot(plotDim.nRow,plotDim.nCol,plotDim.pos(iChan));
                % prepare UserData for current subplot
                spData.chanData = chanData;
                spData.plotPar = plotPar;
                spData.currInd = currInd;
                spData.labels = labels;
                spData.analysis = analysis;
                spData.paths = paths;
                spData.taskType = taskType;
                spData.sigInt.raw = [];
                spData.sigInt.fdr = [];
                % assign UserData of current subplot
                set(erpSpHandle,'UserData',spData);
                set(erpSpHandle,'ButtonDownFcn',{@SubplotCallback,erpSpHandle});
                % plot the current channel
                pHandle = subplotERP('electrode array','channel data',chanData,'plot par',plotPar,'analysis',analysis);
            end
            % prepare place of legend
            h2 = subplot(plotDim.nRow,plotDim.nCol,plotDim.legend); hold on;
            set(subplot(plotDim.nRow,plotDim.nCol,plotDim.legend),'Color',[0.8 0.8 0.8])
            axis off;
            
            % Create Legend
            hL=legend(pHandle,currLabels);
            set(hL,...
                'Position', get(h2,'position'),...
                'Units', get(h2,'Units'),...
                'Box', 'off');
            
            clear h h2 hL;
            
            % Insert UI elements
            UIhandle.vol = uicontrol(fh,'Style', 'checkbox', 'String', 'Voltage',...
                                     'Position', [5 30 80 20]);
            UIhandle.scd = uicontrol(fh,'Style', 'checkbox', 'String', 'SCD',...
                                     'Position', [5 10 80 20]);
            UIhandle.btn = uicontrol(fh,'Style', 'pushbutton', 'String', 'Plot',...
                                     'Position', [90 10 80 40],...
                                     'Callback', {@topoPlot,UIhandle,plotPar.figInd});
            UIhandle.newPlot = uicontrol(fh,'Style', 'pushbutton', 'String', 'New Plot',...
                                         'Position', [175 10 80 40]);
            % Set the Userdata for the button of the current plot (including the data,
            % conditions, plotting parameters etc.)
            topo.plotPar = plotPar;
            topo.analysis = analysis;
            topo.paths = paths.local;
            topo.data = erpAll(:,currInd,:,:);
            topo.conds = plotConds{iCond};
            topo.chanlocs = chanlocs;
            set(UIhandle.btn,'UserData',topo);
            % Set UserData for NewPlot Button
            plotData.erpAll = urerpAll;
            plotData.chanlocs = urchanlocs;
            plotData.plotPar = plotPar;
            plotData.trig = urtrig;
            plotData.analysis = analysis;
            plotData.paths = paths;
            plotData.taskType = taskType;
            plotData.labels = labels;
            set(UIhandle.newPlot,'UserData',plotData);
            set(UIhandle.newPlot,'Callback',{@plotDialog,UIhandle});
        end
    end
    
    if analysis.statsFlag
        %% Plot curves at stats channel with histograms and extract data
        % get position of channels for plotting
        channels2plot = zeros(1,size(plotPar.plotChannelsStat,1));
        for nChannel = 1:size(plotPar.plotChannelsStat,1) 
            for n=1:numel(chanlocs)
                if strcmp(chanlocs(n).labels,plotPar.plotChannelsStat(nChannel,:))
                    channels2plot(nChannel) = n;
                end
            end
        end
        
        % for each specified stats channel
        for iChan = 1:size(channels2plot,2)
            plotPar.currChan = iChan;
            plotPar.currChanLabel = plotPar.plotChannelsStat{iChan};
            for iFig = 1:size(plotPar.plotCondsStat,2)
                fh = figure('WindowButtonDownFcn',@closeFig);
                plotPar.figInd.curr = iFig;
                % determine which cells in plotPar.plotConds belong to the current task
                % and rearrange it
                plotConds = {};
                % get colors
                if plotPar.autoColor
                    % plotColor = distinguishable_colors(max(max(cellfun(@length,plotPar.plotCondsStat)))); 
                    plotColor = distinct_color(max(max(cellfun(@length,plotPar.plotCondsStat)))); 
                else
                    % plotColor = distinguishable_colors(max(trig.color(:,plotPar.plotCondsStatCol(iFig))));
                    plotColor = distinct_color(max(trig.color(:,plotPar.plotCondsStatCol(iFig))));
                end
                for iPlot = 1:size(plotPar.plotCondsStat,1)
                    % go through each cell in plotPar.plotConds
                    if ~isempty(plotPar.plotCondsStat{iPlot,iFig})
                        if any(ismember(taskType,cell2mat(plotPar.plotCondsStat{iPlot,iFig}(end))))
                            plotConds{end+1} = plotPar.plotCondsStat{iPlot,iFig}(1:end-1);
                        end
                    end
                end
                
                % only when there are waves to plot for the current task
                if ~isempty(plotConds)        
                    % get minimal and maximal value for current figure
                    maxVal = getMax('statistics',erpAll,analysis,plotPar,plotConds,labels,chanlocs,iChan,channels2plot);
                    plotPar.yScale(1) = maxVal.erpMin;
                    plotPar.yScale(2) = maxVal.erpMax;
                    if ~isempty(plotPar.compWin)
                        plotPar.yScaleMean(1) = maxVal.meanMin;
                        plotPar.yScaleMean(2) = maxVal.meanMax;
                    else
                        plotPar.yScaleMean = [];
                    end

                    % get subplot dimension parameters
                    plotDim = getDisp('statistics','structure',plotConds,'parameters',plotPar);
                    % initialize result matrix and headers for saving 
                    spData.resAllParent = [];
                    spData.resHeadParent = {};
                    for iPlot = 1:size(plotConds,2)
                        % initialize variables
                        currInd = [];
                        plotPar.currColor = [];
                        plotPar.currStyle = [];
                        % initialize result matrix and header for subplotCallBack saving
                        spData.resAll = [];
                        spData.resHead = {};
                        % get the data to plot 
                        for iCurve = 1:size(plotConds{iPlot},2) 
                            % go through each Curve of the current figure
                            foundLabel = 0;
                            for iIndex = 1:size(labels,1)
                                % check current curve label against each label of the
                                % curves stored in erpAll to get the index of the current
                                % curve label 
                                if strcmp(plotConds{iPlot}(iCurve),labels(iIndex))
                                    foundLabel = 1;
                                    currInd = [currInd iIndex];
                                    if plotPar.autoColor
                                        plotPar.currColor = [plotPar.currColor;  plotColor(iCurve,:);]; 
                                    else
                                        plotPar.currColor = [plotPar.currColor; plotColor(trig.color(iIndex,plotPar.plotCondsStatCol(iFig)),:)];
                                        plotPar.currStyle = [plotPar.currStyle; plotPar.style(iIndex)];
                                    end
                                end
                            end
                            if ~foundLabel
                                % through error if requested curve label is not availabel
                                error([':: Did not find curve label ' plotConds{iPlot}{iCurve}]);
                            end
                        end

                        % get the labels of the current plot
                        currLabels = plotConds{iPlot};
                        chanData = [];            
                        % get the data for the current condition and the current channel
                        for iCurve = 1:size(currInd,2)
                            chanData(:,iCurve) = squeeze(gavr(:,currInd(iCurve),:,channels2plot(iChan)));
                        end
                        % create the subplot for the current channel
                        erpSpHandle = subplot(plotDim.nRow,plotDim.nCol,plotDim.curvePos(iPlot));
                        % prepare UserData for current subplot
                        spData.chanData = chanData;
                        spData.plotPar = plotPar;
                        spData.currInd = currInd;
                        spData.labels = labels;
                        spData.analysis = analysis;
                        spData.paths = paths;
                        spData.taskType = taskType;
                        spData.channelIndex = channels2plot(iChan);

                        if plotPar.runningStat
                            % running statistics (one way RMANOVA) over each time point
                            % in given range for current plot
                            disp([':: Calculating ' plotPar.statTest  ' for plot ' ...
                                  num2str(iPlot) ' in figure ' num2str(iFig)]);
                            sigInt = runningStat(erpAll,spData);
                        else
                            sigInt.raw = [];
                            sigInt.fdr = [];
                        end
                        spData.sigInt = sigInt;
                        % plot the current channel
                        pHandle = subplotERP('statistics','channel data',chanData,'sig',sigInt,'plot par',plotPar,'analysis',analysis);
                        h = legend(pHandle,plotConds{iPlot});
                        set(h, 'Location', 'southwest');
                        
                        if ~isempty(plotPar.compWin)
                            % Plot Bar diagram of components
                            statSpHandle = subplot(plotDim.nRow,plotDim.nCol,plotDim.histPos(iPlot)); 
                            hold on;
                            
                            % for each component window
                            erpMean = zeros(size(plotPar.compWin,1),size(chanData,2));
                            erpErr = erpMean;
                            
                            for iWin = 1:size(plotPar.compWin,1)
                                
                                statWin = [round(((plotPar.compWin(iWin,1)+abs(plotPar.xScale(1)))*analysis.sampRate)/1000) ...
                                           round(((plotPar.compWin(iWin,2)+abs(plotPar.xScale(1)))*analysis.sampRate)/1000)];
                                
                                % get data for each wave
                                for iWave = 1:size(chanData,2)
                                    erpMean(iWin,iWave) = mean(mean(squeeze(erpAll(:,currInd(iWave),statWin(1):statWin(2),channels2plot(iChan))),2));
                                    erpErr(iWin,iWave) = std(mean(squeeze(erpAll(:,currInd(iWave),statWin(1):statWin(2),channels2plot(iChan))),2))/sqrt(size(erpAll,1));
                                    % Save Data in results matrix columns resemble conditions 
                                    spData.resAll = [spData.resAll mean(squeeze(erpAll(:,currInd(iWave),statWin(1):statWin(2),channels2plot(iChan))),2)];
                                    spData.resAllParent = [spData.resAllParent mean(squeeze(erpAll(:,currInd(iWave),statWin(1):statWin(2),channels2plot(iChan))),2)];                                
                                    % write header for current column
                                    spData.resHead = [spData.resHead [char(plotPar.plotCondsStat{iPlot,iFig}(iWave)) 'Win' num2str(iWin)]];
                                    spData.resHeadParent = [spData.resHeadParent [char(plotPar.plotCondsStat{iPlot,iFig}(iWave)) 'Win' num2str(iWin)]];
                                end
                            end 
                            % store erpMean and erpErr values in spData for subplot Callback
                            spData.erpMean = erpMean;
                            spData.erpErr = erpErr;
                            
                            % plot histogram
                            if size(plotPar.comps,1) == 1
                                % plot each bar separately if the there is only one component window
                                % (otherwise there will be problems with the bar
                                % coloring)
                                for iBar = 1:size(erpMean,2)
                                    bar(iBar,erpMean(iBar),'facecolor',plotPar.currColor(iBar,:));
                                end
                            else
                                hB = bar(erpMean);
                            end
                            % color the bars according to the waves in the plots
                            if size(plotPar.comps,1) > 1
                                for iWave = 1:size(erpMean,2)
                                    set(hB(iWave),'facecolor',plotPar.currColor(iWave,:));
                                end
                            end

                            numgroups = size(erpMean, 1);
                            numbars = size(erpMean, 2);
                            
                            groupwidth = min(0.8, numbars/(numbars+1.5));

                            % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
                            xVec = [];
                            for iBar = 1:numbars
                                if size(plotPar.comps,1) == 1
                                    x = iBar;
                                    xVec = [xVec x];
                                else                        
                                    x = (1:numgroups) - groupwidth/2 + (2*iBar-1) * groupwidth / (2*numbars); % Aligning error bar with individual bar
                                    xVec = [xVec x];
                                end
                                errorbar(x, erpMean(:,iBar), erpErr(:,iBar), 'k', 'linestyle', 'none');
                            end
                            % reshape the bar positions according to the component windows
                            if size(erpMean,1) > 1
                                xVec = reshape(xVec,size(erpMean,1),size(erpMean,2));
                            else
                                set(gca,'XLim',[(xVec(1) - 1) (xVec(end) + 1)]);
                            end
                            % get the central position (centarl bar) of each component windows
                            xTickPos = zeros(1,size(xVec,1));
                            for iWin = 1:size(xVec,1)
                                xTickPos(1,iWin) = median(xVec(iWin,:));
                            end
                            set(gca,'XTick',xTickPos);
                            set(gca,'XTickLabel',plotPar.winNames);
                            set(gca,'YLim',[(maxVal.meanMin - plotPar.yOverhead)  (maxVal.meanMax + plotPar.yOverhead)]);
                            set(get(gca,'YLabel'),'String','Voltage (micro Volts)');
                            set(statSpHandle,'ButtonDownFcn',{@SubplotCallback,erpSpHandle});
                        
                            % add savebutton
                            UIsave = uicontrol(fh,'Style', 'pushbutton', 'String', 'Save Data',...
                                               'Position', [5 10 100 40]);
                            
                            % Assigning relevant data for UI UserData
                            set(UIsave,'UserData',spData);
                            set(UIsave,'Callback',{@saveRes,UIsave});
                        end
                        % assign UserData of current subplot
                        set(erpSpHandle,'UserData',spData);
                        set(erpSpHandle,'ButtonDownFcn',{@SubplotCallback,erpSpHandle});
                    end % for iPlot
                    
                    if isempty(plotPar.compWin)
                        UIhandle.newPlot = uicontrol(fh,'Style', 'pushbutton', 'String', 'New Plot',...
                                                     'Position', [5 10 100 40]);
                    else
                        UIhandle.newPlot = uicontrol(fh,'Style', 'pushbutton', 'String', 'New Plot',...
                                                     'Position', [110 10 80 40]);
                    end
                    % Set UserData for NewPlot Button
                    plotData.erpAll = urerpAll;
                    plotData.chanlocs = urchanlocs;
                    plotData.plotPar = plotPar;
                    plotData.trig = urtrig;
                    plotData.analysis = analysis;
                    plotData.paths = paths;
                    plotData.taskType = taskType;
                    plotData.labels = labels;
                    set(UIhandle.newPlot,'UserData',plotData);
                    set(UIhandle.newPlot,'Callback',{@plotDialog,UIhandle});
                else
                    % if plotConds is empty (because there are no waves for
                    % the current task - close the current figure
                    close(fh);
                end
            end % for iFig
        end % for iChan
    end % if plotPar.statsFlag
end

%% CALLBACK FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% execute function when mouseclick on subplot
function SubplotCallback(src,eventdata,erpSpHandle,fh)

    spData = get(erpSpHandle,'UserData');
    % evaluate display method
    switch spData.plotPar.singleDispMode
      case 1
        % open one more figure then scheduled by plot_erp and override this every
        % time SubplotCallback is called
        singleFig = figure(1000);
        set(singleFig,'name','Zoomed Figure');
        clf;
        set(gcf,'ButtonDownFcn',@closeFig,'Color',[0.8 0.8 0.8]);              
      case 2
        % open a new figure for each click
        figure
        set(gcf,'ButtonDownFcn',@closeFig,'Color',[0.8 0.8 0.8]);
      otherwise
        error([':: Invalid option for plotPar.singleDispMode: ' num2str(plotPar.singleDispMode)]);
    end
    
    % plot the clicked subplot in this new figure with statistics
    singleERP(spData);
end

function topoPlot(hObject,eventdata,UIhandle,figInd)
    
% get the User data
    topo = get(hObject,'UserData');
    % get states of checkboxes
    topo.vol = get(UIhandle.vol,'value');
    topo.scd = get(UIhandle.scd,'value');
    % promt list when clicking plot to choose waves 
    [selection,ok]=listdlg('ListString',topo.conds,'PromptString','Choose waves.');
    % select the waves
    topo.conds = topo.conds(selection);
    % select the data
    topo.data = topo.data(:,selection,:,:);
    % plot topographies of selected data
    if topo.vol | topo.scd
        preptopo(topo);
    else
        errordlg('Choose at least one plotting method.');
    end
end

function plotDialog(src,evnt,UIhandle)
% get the User data
    plotData = get(UIhandle.newPlot,'UserData');
    % clear plotConds and plotCondsStat
    figNum = 0;
    plotData.plotPar.plotConds = {};
    plotData.plotPar.plotCondsStat = {};
    % set defaults for Plotting Maps
    plotData.analysis.plotERPflag = 0;
    plotData.analysis.statsFlag = 0;
    % Initialize plotstring
    plotString = [];
    plotStringSep = '----------------------------------------';
    
    UIhandle.menuH = dialog('WindowStyle', 'normal');
    % TextBox for scheduled fugures
    UIhandle.txtBox = uicontrol('Parent',UIhandle.menuH,...
                                'Units','normalized',...
                                'Position',[0.2,0.2,0.75,0.75],...
                                'BackgroundColor', [1 1 1],...
                                'Style','edit',...
                                'Max',100,...
                                'Enable','inactive',...
                                'String',plotString);
    UIhandle.array = uicontrol(UIhandle.menuH,'Style', 'checkbox', 'String', 'Electrode Array',...
                               'Position', [5 50 120 20]);
    UIhandle.stats = uicontrol(UIhandle.menuH,'Style', 'checkbox', 'String', 'Statistics',...
                               'Position', [5 70 100 20]);
    UIhandle.autoColor = uicontrol(UIhandle.menuH,'Style', 'checkbox', 'String', 'Auto Color',...
                               'Position', [195 5 100 20]);
    
    uicontrol(UIhandle.menuH,'Style', 'pushbutton', 'String', 'Close',...
              'Position', [5,396,50,20],...
              'Callback', {@closeDialog,UIhandle.menuH});
    uicontrol(UIhandle.menuH,'Style', 'pushbutton', 'String', 'Plot',...
              'Position', [5 5 80 40],...
              'BackgroundColor', [0.8 0.2 0],...
              'Callback', {@newPlot,UIhandle});
    uicontrol(UIhandle.menuH,'Style', 'pushbutton', 'String', 'Add Waves',...
              'Position', [90 5 100 40],...
              'Callback', {@addWaves,UIhandle});
    

    
    function addWaves(src,evnt,UIhandle)
    % promt list when clicking plot to choose waves 
        if get(UIhandle.array,'value') | get(UIhandle.stats,'value')
            [selection,ok]=listdlg('ListString',plotData.labels,...
                                   'PromptString','Choose waves.',...
                                   'ListSize', [300 300]);
        else
            errordlg('Choose at least one plotting method.');
        end
        % evaluate plotting methods
        if get(UIhandle.array,'value') & get(UIhandle.stats,'value') & ~isempty(selection)
            figNum = figNum + 1;
            plotData.analysis.plotERPflag = 1;
            plotData.analysis.statsFlag = 1;
            plotHeader = [plotStringSep sprintf('\n')...
                          'Figure ' num2str(figNum)...
                          ' (Array and Statistics)' sprintf('\n')];
            % Add Selected curves to plotConds
            plotData.plotPar.plotConds{end+1} = plotData.labels(selection)';
            plotData.plotPar.plotCondsStat{end+1} = plotData.labels(selection)';
            % Add TextBox to show the scheduled plot
            plotString = [plotString plotHeader plotData.plotPar.plotConds{end} plotStringSep];
            set(UIhandle.txtBox,'String',plotString);
            % Add the task type
            plotData.plotPar.plotConds{end} = [plotData.plotPar.plotConds{end} plotData.taskType];
            plotData.plotPar.plotCondsStat{end} = [plotData.plotPar.plotCondsStat{end} plotData.taskType];
        elseif get(UIhandle.array,'value') & ~isempty(selection)
            figNum = figNum + 1;
            plotData.analysis.plotERPflag = 1;
            plotHeader = [plotStringSep sprintf('\n')...
                          'Figure ' num2str(figNum)...
                          ' (Electrode Array)' sprintf('\n')];
            % Add Selected curves to plotConds
            plotData.plotPar.plotConds{end+1} = plotData.labels(selection)';
            % Add TextBox to show the scheduled plot
            plotString = [plotString plotHeader plotData.plotPar.plotConds{end} plotStringSep];
            set(UIhandle.txtBox,'String',plotString);
            % add the task type
            plotData.plotPar.plotConds{end} = [plotData.plotPar.plotConds{end} plotData.taskType];
        elseif get(UIhandle.stats,'value') & ~isempty(selection)
            figNum = figNum + 1;
            plotData.analysis.statsFlag = 1;
            plotHeader = [plotStringSep sprintf('\n')...
                          'Figure ' num2str(figNum)...
                          ' (Statistics)' sprintf('\n')];
            % Add Selected curves to plotConds
            plotData.plotPar.plotCondsStat{end+1} = plotData.labels(selection)';
            % Add TextBox to show the scheduled plot
            plotString = [plotString plotHeader plotData.plotPar.plotCondsStat{end} plotStringSep];
            set(UIhandle.txtBox,'String',plotString);
            % add the task type
            plotData.plotPar.plotCondsStat{end} = [plotData.plotPar.plotCondsStat{end} plotData.taskType];
        else

            return;
        end
        set(UIhandle.newPlot,'UserData',plotData);
    end
end

function newPlot(src,evnt,UIhandle)
% retrieve plotting Data
    plotData = get(UIhandle.newPlot,'UserData');
    % Assign auto color value
    plotData.plotPar.autoColor = get(UIhandle.autoColor,'value');
    
    % call plot functions
    plot_erp(plotData.erpAll,plotData.chanlocs,plotData.plotPar,plotData.trig,...
             plotData.analysis,plotData.paths,plotData.taskType);
    % close Dialog box
    delete(UIhandle.menuH);
end

function saveRes(src,evnt,UIsave)
    
% retrieve UserData
    spData = get(UIsave,'UserData');
    % open save Dialog
    [fileName,pathName] = uiputfile([spData.paths.resDir spData.paths.taskLabel{spData.taskType} '.txt']);
    
% save the data
fid = [pathName fileName];
% the engine
txt=sprintf('%s\t',spData.resHeadParent{:});
txt(end)='';
disp(':: Saving Data.');
dlmwrite(fid,txt,'');
dlmwrite(fid,spData.resAllParent,'-append','delimiter','\t');

end


function closeDialog(src,evnt,hObject)
   delete(hObject) 
end

function closeFig(src,evnt)
% Press Ctrl Click to close
    if strcmp(get(gcf,'SelectionType'),'alt')
    close(gcf); 
    else
        % disp(':: Ctrl-Click to close');
    end
end

