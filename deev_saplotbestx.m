function stop = deev_saplotbestx(options,optimvalues,flag)
%SAPLOTBESTX PlotFcn to plot best X value.
%   STOP = SAPLOTBESTX(OPTIONS,OPTIMVALUES,FLAG) where OPTIMVALUES is a
%   structure with the following fields:
%              x: current point 
%           fval: function value at x
%          bestx: best point found so far
%       bestfval: function value at bestx
%    temperature: current temperature
%      iteration: current iteration
%      funccount: number of function evaluations
%             t0: start time
%              k: annealing parameter
%
%   OPTIONS: The options structure created by using SAOPTIMSET
%
%   FLAG: Current state in which PlotFcn is called. Possible values are:
%           init: initialization state
%           iter: iteration state
%           done: final state
%
%   STOP: A boolean to stop the algorithm.
%
%   Example:
%    Create an options structure that will use SAPLOTBESTX as the plot
%    function
%     options = saoptimset('PlotFcns',@saplotbestx);

%   Copyright 2006-2010 The MathWorks, Inc.
ub = evalin('base','ub');
stop = false;
switch flag
    case 'init'
        set(gca,'xlimmode','manual','zlimmode','manual', ...
            'alimmode','manual')
        substr = evalin('base','substr');
        evtStr = evalin('base','evtStr');
        title(sprintf('Best point-%s-%s',substr,evtStr),'interp','none')
        Xlength = numel(optimvalues.bestx);
        xlabel(sprintf('Number of variables (%i)',Xlength),'interp','none');
        ylabel('Best point, pct of max','interp','none');
        plotBestX = bar(optimvalues.bestx(:)./ub');
        set(plotBestX,'Tag','saplotbestx');
        set(plotBestX,'edgecolor','none')
        set(gca,'xlim',[0,1 + Xlength])
        hold on 
        for it = 1:Xlength
            text(plotBestX.XData(it)-.25, .5, num2str(optimvalues.bestx(it),'%.02f'), 'color', [1 0 0],'Tag',['x' num2str(it)]);
        end
    case 'iter'
        plotBestX = findobj(get(gca,'Children'),'Tag','saplotbestx');
        set(plotBestX,'Ydata',optimvalues.bestx(:)./ub')
        x = nan(1,4);
        x(1) = optimvalues.bestx(1); x(3) = optimvalues.bestx(3);
        x(2) = x(1)+optimvalues.bestx(2); x(4) = x(3)+optimvalues.bestx(4);
        Xlength = numel(optimvalues.bestx);
        str = {'s','e','l','u'};
        for ih = 1:Xlength
            h = findobj(gcf,'Tag',['x' num2str(ih)]);
            set(h,'String',[str{ih} ':' num2str(x(ih),'%.02f')]); 
        end 
        ylim([0 1]);
end