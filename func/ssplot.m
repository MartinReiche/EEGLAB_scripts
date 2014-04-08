% This function takes subjects X trigger X sample data and plots
% one separate line for each subject in one plot per trigger type to
% observe single subject ERPs per trigger
%
% Copyright (c) 2014 Martin Reiche, Carl-von-Ossietzky-University Oldenburg
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

function ssplot(spData)

% switch auto dimensions for subplot structure on
    spData.plotPar.autoDim = 1;
    dispPars = getDisp('electrode array','parameters',spData.plotPar,'structure',[1:size(spData.chanSingleData,2)]);

    % open new figure
    fh = figure;
    set(fh,'ButtonDownFcn',@closeFig,'Color',[0.8 0.8 0.8]);

    % get color values 
    spData.plotPar.currColor = distinct_color(size(spData.chanSingleData,1));
    

    
    for iPlot = 1:size(spData.chanSingleData,2)
        % create one subplot for each trigger
        erpSpHandle = subplot(dispPars.nRow,dispPars.nCol,dispPars.pos(iPlot));
        % select data of current condition
        chanData = squeeze(spData.chanSingleData(:,iPlot,:))';
        % plot all subject ERPs of current condition in one plot
        subPlotHandle = subplotERP('single subject','channel data',chanData,'plot par',spData.plotPar,'analysis',spData.analysis,'sig',spData.sigInt);
        spData.chanData = chanData;
        tH = title(spData.labels(spData.currInd(iPlot)));
        spData.title = spData.labels(spData.currInd(iPlot));
        
        YLim = get(erpSpHandle,'YLim');
        tLoc = get(tH,'Position');
        tLoc(2) = YLim(1) - 0.3;
        set(tH,'Position',tLoc);
        
        set(erpSpHandle,'UserData',spData);
        set(erpSpHandle,'ButtonDownFcn',{@SubplotCallback,erpSpHandle});
    end

    % prepare legend labels
    legLabels = {};
    for iSub = 1:numel(spData.analysis.subjects)
        legLabels{iSub} = ['Subject ' num2str(spData.analysis.subjects(iSub))];
    end
    
    % prepare legend position
    h2 = subplot(dispPars.nRow,dispPars.nCol,dispPars.legend); hold on;
    set(subplot(dispPars.nRow,dispPars.nCol,dispPars.legend),'Color',[0.8 0.8 0.8])
    axis off;
    % Create Legend
    hL=legend(subPlotHandle,legLabels);
    set(hL,...
        'Position', get(h2,'position'),...
        'Units', get(h2,'Units'),...
        'Box', 'off');
    
    set(erpSpHandle,'UserData',spData);
    set(erpSpHandle,'ButtonDownFcn',{@SubplotCallback,erpSpHandle});

    
end

%% CALLBACK FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% execute function when mouseclick on subplot
function SubplotCallback(src,eventdata,erpSpHandle,fh)
    
    spData = get(erpSpHandle,'UserData');
    spData.type = 'singleSub';
    spData.plotPar.singleScaleAuto = 1;
    figure
    set(gcf,'ButtonDownFcn',@closeFig,'Color',[0.8 0.8 0.8]);

    spData.plotPar.compWin = [];
    % plot the clicked subplot in this new figure with statistics
    singleERP(spData);
end

function closeFig(src,evnt)
% Press Ctrl Click to close
    if strcmp(get(gcf,'SelectionType'),'alt')
        close(gcf); 
    else
        % disp(':: Ctrl-Click to close');
    end
end
