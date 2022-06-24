classdef basePlot < handle
    properties
        calcObj
        lineWidth   = 2.0
        colorsFull  = colororder
        marksOrder  = {'+','o','*','x','s','d','^','v','>','<','p','h'}
        marksLatex  = {'+','$\circle$','*','x','$\square$','$\diamond$',...
            '$\triangle$','$\triangledown$','$\triangleright$',...
            '$\triangleleft$','$\star$','$\hexagram$'}
        marksLW     = 1.0
        marksSize   = 6.0
    end
    properties (SetAccess=protected)
        figs    matlab.ui.Figure
        axs     matlab.graphics.axis.Axes
    end
    properties (Dependent)
        nFigs
        nAxs
        nCol
        nMarks
    end
    methods
        function obj = basePlot(calcObj)
            obj.calcObj = calcObj;
        end
        function val = get.nMarks(obj)
            val = length(obj.marksOrder);
        end
        function val = get.nCol(obj)
            val = size(obj.colorsFull,1);
        end
        function val = get.nFigs(obj)
            val = length(obj.figs);
        end
        function val = get.nAxs(obj)
            val = length(obj.axs);
        end
        function [fg,ax]=prePlot(~)
            % PREPLOT Basic plot initializations

            %% Create basic figure object
            fg=figure;ax=axes(fg);
            hold(ax,'on');
            ax.Color=[1.9*0.5 1.9*0.5 1.9*0.5];
        end

        function postPlot(obj,fg,ax,Args)
            arguments
                obj
                fg
                ax
                Args.File
            end
            % POSTPLOT Finish plot

            hold(ax,'off');
            %% Check if figure already included
            found = false;
            for curr_fg = obj.figs
                if fg == curr_fg
                    found = true;
                    break
                end
            end
            %% Add to tracking list
            if ~found
                obj.figs(obj.nFigs + 1) = fg;
                obj.axs(obj.nAxs + 1) = ax;
            end
            %% Export to file
            if isfield(Args,'File')
                title = ax.Title.String;
                ax.Title.String=[];
                subArgs=cell(1,0);
                [path,name,ext] = fileparts(Args.File);
                switch ext
                    case '.png'
                        subArgs = [subArgs(:)',{'Resolution', 600}];
                    case '.pdf'
                        subArgs = [subArgs(:)',{'ContentType','vector', 'BackgroundColor','none'}];
                end
                exportgraphics(fg,Args.File,subArgs{:});
                ax.Title.String=title;
                savefig(fg,[path name '.fig'],'compact');
            end
        end

        function [x,y,meta]=breakLines(~,x,y,Args)
            arguments
                ~
                x
                y
                Args.meta
                Args.yTresh = 1
            end
            % BREAKLINES Handle discontinuities

            ix = 2;
            if isfield(Args,'meta')
                meta = Args.meta;
            end
            while ix ~= length(x)
                ty = squeeze(y(ix,:,:));
                ty_prev = squeeze(y(ix-1,:,:));
                if isfield(Args,'meta')
                    tmeta = squeeze(meta(ix,:,:));
                    tmeta_prev = squeeze(meta(ix-1,:,:));
                end
                flag = abs(ty - ty_prev) > Args.yTresh;
                if sum(flag(:))
                    x_mid = (x(ix-1) + x(ix)) / 2;
                    y_mid = (ty + ty_prev) / 2;
                    y_mid(flag(:)) = nan;
                    x = [x(1:ix-1) x_mid x(ix:end)];
                    y((ix:end)+1,:,:) = squeeze(y(ix:end,:,:));
                    y(ix,:,:) = y_mid;
                    if isfield(Args,'meta')
                        meta_mid = (tmeta + tmeta_prev) / 2;
                        meta_mid(flag(:)) = nan;
                        meta((ix:end)+1,:,:) = squeeze(meta(ix:end,:,:));
                        meta(ix,:,:) = meta_mid;
                    end
                    ix = ix + 2;
                else
                    ix = ix + 1;
                end
            end
        end
        function ax2 = zoomin(~,ax,area,position,Args)
            arguments
                ~
                ax
                area        Plot.rect
                position    Plot.rect
                Args.XTick  = []
                Args.YTick  = []
                Args.Corners    = {'top left','bot right'}
            end
            % ZOOMIN add a zoomed in plot on area
            
            fig = ax.Parent;
            ax2 = copyobj(ax,fig);
            ax2.Position = double(position);
            ax2.XLim = area.x;
            ax2.YLim = area.y;
            ax2.XTick = Args.XTick;
            ax2.YTick = Args.YTick;
            ax2.XLabel = [];
            ax2.YLabel = [];
            rectangle(ax,'Position',double(area))
            axPos = area.normalize(ax);
            for iC=1:length(Args.Corners)
                switch Args.Corners{iC}
                    case 'top left'
                        annotation(fig,'line',[axPos.x(1) position.x(1)],...
                            [axPos.y(1) position.y(2)],'Color','k');
                    case 'top right'
                        annotation(fig,'line',[axPos.x(2) position.x(2)],...
                            [axPos.y(1) position.y(2)],'Color','k');
                    case 'bot left'
                        annotation(fig,'line',[axPos.x(1) position.x(1)],...
                            [axPos.y(2) position.y(1)],'Color','k');
                    case 'bot right'
                        annotation(fig,'line',[axPos.x(2) position.x(2)],...
                            [axPos.y(2) position.y(1)],'Color','k');
                    otherwise
                        error('Unknown corner');
                end
            end
        end
    end
end