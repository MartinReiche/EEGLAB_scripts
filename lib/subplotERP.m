% Takes a Matrix of ERP Data (line represent sampling points, columns
% represent curves). Depending on the columns, plots the curves in one plot
% with predefined colours, adds labels and adjusts the graphics.
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

function plotHandle = subplotERP(method,varargin)

%% Evaluate input
for iArg = 1:2:length(varargin)
    switch lower(varargin{iArg})
      case 'channel data'
        chanData = varargin{iArg+1};
      case 'sig'
        sigInt = varargin{iArg+1};
      case 'plot par'
        plotPar = varargin{iArg+1};
      case 'analysis'
        analysis = varargin{iArg+1};
      otherwise
        error([':: ' varargin{iArg} ' is not a valid option']);
    end
end


%% get x scaling
% get baseline in Seconds 
baselineMs = abs(analysis.erpWin(1));
% get time resolution
timeRes = 1000/analysis.sampRate;
% number of samples
noSamp = size(chanData,1);
% get number of samples in the baseline
baselinePnts = ceil(baselineMs/timeRes);
% create vector of sample points in seconds
x = [-baselinePnts*timeRes:timeRes:noSamp*timeRes-baselinePnts*timeRes-timeRes];

% evaluate method
switch lower(method)
    %% plot electrode array for given data
  case 'electrode array'
    %% PLOTTING
    hold on;
    % add component boxes
    for nComp = 1:size(plotPar.comps,1)
        if plotPar.compWin(nComp,2)-plotPar.compWin(nComp,1) > 0
            % add marking of current component (window)
            rectangle('Position',[plotPar.compWin(nComp,1) plotPar.yScale(1)...
                                plotPar.compWin(nComp,2)-plotPar.compWin(nComp,1)...
                                plotPar.yScale(2)-plotPar.yScale(1)],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
            % add component name of current window
            text(plotPar.compWin(nComp,1),plotPar.yScale(1)-0.25*plotPar.yCoef,plotPar.comps{nComp},'FontSize',11);

        end
    end
    
    % Draw the baseline interval
    if plotPar.drawBaseLine && analysis.rmBase
        % mark baseline
        rectangle('Position',[plotPar.baseWin(1) plotPar.yScale(1)...
                            plotPar.baseWin(2)-plotPar.baseWin(1)...
                            plotPar.yScale(2)-plotPar.yScale(1)],'FaceColor',[0.9 ...
                            0.9 0.9],'EdgeColor','none');
        % plot baseline label
        text(plotPar.baseWin(1)+10,plotPar.yScale(1)-0.15*plotPar.yCoef,'Baseline','FontSize',11);
    end

    % add name of the electrode to the current plot
    text(plotPar.xScale(1)+10,plotPar.yScale(1)+0.45*plotPar.yCoef,plotPar.plotChannels{plotPar.currChan},'FontSize',16);
    % add labels
    xlabel('Latency (ms)');
    ylabel(['Amplitude (microVolts)']);
    % adjust Axes
    axis([plotPar.xScale(1) plotPar.xScale(2) plotPar.yScale(1) plotPar.yScale(2)]);
    % set the correct X and Y ticks and reverse ordinate
    % get the XTicks
    secondTick = ((ceil(plotPar.xScale(1)/plotPar.xCoef)-((plotPar.xScale(1)/plotPar.xCoef)))*plotPar.xCoef)+plotPar.xScale(1);

    % construct the x scale
    xRange = [secondTick:plotPar.xCoef:plotPar.xScale(2)];
    % check whether xRange adds up to the end of the spectra
    if xRange(end) ~= plotPar.xScale(2) 
        % if not, add it
        xRange = [xRange plotPar.xScale(2)];
    end
    if xRange(1) ~= plotPar.xScale(1) 
        % if not, add it
        xRange = [plotPar.xScale(1) xRange];
    end
     
    set(gca,'Xtick',xRange,'Ytick',plotPar.yScale(1):plotPar.yCoef:plotPar.yScale(2),...
            'YDir','reverse','box','off','FontSize',11);
    % add zero lines
    set(line(plotPar.xScale,[0 0]),'Color',[0 0 0],'LineStyle','-','linewidth',0.5);
    set(line([0 0],plotPar.yScale),'Color',[0 0 0],'LineStyle','-','linewidth',0.5);
    % add grid
    if plotPar.grid
        grid on;
    end
    % start plotting
    for iCurve = 1:size(chanData,2)
        % plot current curve
        plotHandle(iCurve) = plot(x,chanData(:,iCurve),'color',plotPar.currColor(iCurve,:),...
                                  'LineWidth',plotPar.lineWidth,'LineStyle',plotPar.currStyle{iCurve});
    end
    hold off;
    %% plot electrode array for given data
  case 'statistics'
    %% PLOTTING
    hold on;
    
    % Draw the baseline interval
    if plotPar.drawBaseLine && analysis.rmBase
        % mark baseline
        rectangle('Position',[plotPar.baseWin(1) plotPar.yScale(1)...
                            plotPar.baseWin(2)-plotPar.baseWin(1)...
                            plotPar.yScale(2)-plotPar.yScale(1)],'FaceColor',[0.9 ...
                            0.9 0.9],'EdgeColor','none');
        % plot baseline label
        text(plotPar.baseWin(1)+10,plotPar.yScale(1)-0.15*plotPar.yCoef,'Baseline','FontSize',11);
    end

    % add component boxes
    for nComp = 1:size(plotPar.comps,1)
        if plotPar.compWin(nComp,2)-plotPar.compWin(nComp,1) > 0
            % add marking of current component (window)
            rectangle('Position',[plotPar.compWin(nComp,1) plotPar.yScale(1)...
                                plotPar.compWin(nComp,2)-plotPar.compWin(nComp,1)...
                                plotPar.yScale(2)-plotPar.yScale(1)],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
            % add component name of current window
            text(plotPar.compWin(nComp,1),plotPar.yScale(1)-0.15*plotPar.yCoef,plotPar.comps{nComp},'FontSize',11);
            
        end
    end
    
    
    % draw significant intervals
    if plotPar.runningStat && ~isempty(sigInt)
        % find start and end of significant intervals
        foundStart = 0;
        intCount = 0;
        diffSig = [diff(sigInt) == 1 0];
        for iInt = 1:size(diffSig,2)
            if diffSig(1,iInt) && ~foundStart
                % increase the interval counter
                intCount = intCount + 1;
                % save start of current interval
                intBound(intCount,1) = sigInt(iInt);
                foundStart = 1;
            elseif ~diffSig(1,iInt) && foundStart
               % save end of current interval
               intBound(intCount,2) = sigInt(iInt);
               foundStart = 0;
            end
        end
        % draw significant intervalls
        for iInt = 1:size(intBound,1)

            % for each intervall
            % determine start end end of current intervall in ms relative to epoch start
            intStart = intBound(iInt,1) * timeRes + analysis.erpWin(1);
            intDur = (intBound(iInt,2) * timeRes +  analysis.erpWin(1))-intStart;
            % draw current intervall
            rectangle('Position',[intStart -plotPar.yScale(2)*0.05 intDur plotPar.yScale(2)*0.05],...
                      'FaceColor',[1 0.5 0.5],'EdgeColor','none');
        end
    end


    % draw tones 
    if plotPar.drawStim 
        % check for drawn stimuli latencies and epoch boundaries
        if (plotPar.stim(1,1) < plotPar.xScale(1)) || (plotPar.stim(end,1) > plotPar.xScale(2))
            warning('Latencies for drawn stimuli exceed epoch boundaries. Skipping.')
        else
            % for each drawn stimulus 
            for iTone = 1:size(plotPar.stim,1)
                
                rectangle('Position',[plotPar.stim(iTone,1) plotPar.stim(iTone,2) plotPar.stimDur 0.2],...
                                    'FaceColor',[0.5 0.5 0.5],'EdgeColor','none');
                line([plotPar.stim(iTone,1) plotPar.stim(iTone,1)],[plotPar.yScale(1) plotPar.yScale(2)],...
                     'LineStyle','--','Color',[0.5 0.5 0.5])
            end
        end
    end

    % add name of the electrode to the current plot
    text(plotPar.xScale(1)+10,plotPar.yScale(1)+0.25*plotPar.yCoef,plotPar.plotChannelsStat{plotPar.currChan},'FontSize',16);
    % add labels
    xlabel('Latency (ms)');
    ylabel(['Amplitude (microVolts)']);
    % adjust Axes
    axis([plotPar.xScale(1) plotPar.xScale(2) plotPar.yScale(1) plotPar.yScale(2)]);
    secondTick = ((ceil(plotPar.xScale(1)/plotPar.xCoef)-((plotPar.xScale(1)/plotPar.xCoef)))*plotPar.xCoef)+plotPar.xScale(1);

    % construct the x scale
    xRange = [secondTick:plotPar.xCoef:plotPar.xScale(2)];
    % check whether xRange adds up to the end of the spectra
    if xRange(end) ~= plotPar.xScale(2) 
        % if not, add it
        xRange = [xRange plotPar.xScale(2)];
    end
    if xRange(1) ~= plotPar.xScale(1) 
        % if not, add it
        xRange = [plotPar.xScale(1) xRange];
    end
        
    set(gca,'Xtick',xRange,'Ytick',plotPar.yScale(1):plotPar.yCoef:plotPar.yScale(2),...
            'YDir','reverse','box','off','FontSize',11);
    % add zero lines
    set(line(plotPar.xScale,[0 0]),'Color',[0 0 0],'LineStyle','-','linewidth',0.5);
    set(line([0 0],plotPar.yScale),'Color',[0 0 0],'LineStyle','-','linewidth',0.5);
    % add grid
    if plotPar.grid
        grid on;
    end
    % start plotting
    for iCurve = 1:size(chanData,2)
        % plot current curve
        plotHandle(iCurve) = plot(x,chanData(:,iCurve),'color',plotPar.currColor(iCurve,:),...
                                  'LineWidth',plotPar.lineWidth,'LineStyle',plotPar.currStyle{iCurve});
    end
    hold off;
      %% Plot single channel ERP with stats when clicked
  case 'single'
    %% PLOTTING
    hold on;
    
    % determine Y Axis scaling
    if plotPar.singleScaleAuto 
        yScale(1) = floor(min(min(chanData)));
        yScale(2) = ceil(max(max(chanData)));
    else
        yScale(1) = plotPar.yScale(1);
        yScale(2) = plotPar.yScale(2);
    end

    % add component boxes
    for nComp = 1:size(plotPar.comps,1)
        if plotPar.compWin(nComp,2)-plotPar.compWin(nComp,1) > 0
            % add marking of current component (window)
            rectangle('Position',[plotPar.compWin(nComp,1) yScale(1)...
                                plotPar.compWin(nComp,2)-plotPar.compWin(nComp,1)...
                                yScale(2)-yScale(1)],'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
            % add component name of current window
            text(plotPar.compWin(nComp,1),yScale(1)-0.15*plotPar.yCoef,plotPar.comps{nComp},'FontSize',11);

        end
    end
    
    if plotPar.drawBaseLine && analysis.rmBase
        % mark baseline
        rectangle('Position',[plotPar.baseWin(1) plotPar.yScale(1)...
                            plotPar.baseWin(2)-plotPar.baseWin(1)...
                            plotPar.yScale(2)-plotPar.yScale(1)],'FaceColor',[0.9 ...
                            0.9 0.9],'EdgeColor','none');
        % plot baseline label
        text(plotPar.baseWin(1)+10,plotPar.yScale(1)-0.15*plotPar.yCoef,'Baseline','FontSize',11);
    end


    % add name of the electrode to the current plot
    text(plotPar.xScale(1)+10,yScale(1)+0.45*plotPar.yCoef,plotPar.currChanLabel,'FontSize',16);
    % add labels
    xlabel('Latency (ms)');
    ylabel(['Amplitude (microVolts)']);
    % adjust Axes
    axis([plotPar.xScale(1) plotPar.xScale(2) yScale(1) yScale(2)]);
    % set the correct X and Y ticks and reverse ordinate
    secondTick = ((ceil(plotPar.xScale(1)/plotPar.xCoef)-((plotPar.xScale(1)/plotPar.xCoef)))*plotPar.xCoef)+plotPar.xScale(1);
    
    % construct the x scale
    xRange = [secondTick:plotPar.xCoef:plotPar.xScale(2)];
    % check whether xRange adds up to the end of the spectra
    if xRange(end) ~= plotPar.xScale(2) 
        % if not, add it
        xRange = [xRange plotPar.xScale(2)];
    end
    if xRange(1) ~= plotPar.xScale(1) 
        % if not, add it
        xRange = [plotPar.xScale(1) xRange];
    end
    
    set(gca,'Xtick',xRange,'Ytick',plotPar.yScale(1):plotPar.yCoef:plotPar.yScale(2),...
            'YDir','reverse','box','off','FontSize',11);
    set(line(plotPar.xScale,[0 0]),'Color',[0 0 0],'LineStyle','-','linewidth',0.5);
    set(line([0 0],yScale),'Color',[0 0 0],'LineStyle','-','linewidth',0.5);
    % add grid
    if plotPar.grid
        grid on;
    end
    % start plotting
    for iCurve = 1:size(chanData,2)
        % plot current curve
        plotHandle(iCurve) = plot(x,chanData(:,iCurve),'color',plotPar.currColor(iCurve,:),...
                                  'LineWidth',plotPar.lineWidth,'LineStyle',plotPar.currStyle{iCurve});
    end
    hold off;
  otherwise
    error([':: ' method ' is not a valid option for method']);
end
end
