classdef TwoLevel < Plot.baseFloquet
    properties (SetAccess=protected)
        colors2     = {'red','blue'}
    end
    methods
        function obj=TwoLevel(calcObj)
            obj@Plot.baseFloquet(calcObj);
        end
%         function PlotAll(obj,Args)
%             arguments
%                 obj
%                 Args.Exact
%                 Args.FHF
%                 Args.NSCF
%             end
%             % PLOTALL Plot all default
% 
%             
%         end
        function [fg,ax] = Variational_th(obj,Res,Args,Args2)
            arguments
                obj
                Res
                Args.filter_nconv   logical =true
                Args.PlotEnpoint    logical =true
                Args2.File
            end
            % VARIATIONAL_OVERLAP Plot the 

            %% Initialize basic figures
            [fg,ax]=obj.prePlot;
            %% Plot Data
            iPlot = 0;
            steps = 1:length(Res(1).steps);
            steps = steps - 1;
            max_steps = 0;
            for iN=1:length(Res)
                if Args.filter_nconv && ~Res(iN).conv; continue; end
                iPlot = iPlot + 1;
                iCol = mod(iPlot - 1, obj.nCol) + 1;
                iMark = mod(iPlot - 1, obj.nMarks) + 1;
                pData = [Res(iN).steps.th];
                ind = find(isnan(pData),1,"first");
                if isempty(ind); ind = steps(end) + 1; end
                max_steps = max(max_steps,ind-1);
                plot(ax,steps,pData,...
                    LineWidth=obj.marksLW,MarkerSize=obj.marksSize,...
                    LineStyle='-',Color=obj.colorsFull(iCol,:),...
                    Marker=obj.marksOrder{iMark});
            end
            % Replot endpoint
            if Args.PlotEnpoint
                for iN=1:length(Res)
                    if Args.filter_nconv && ~Res(iN).conv; continue; end
                    iPlot = iPlot + 1;
                    iCol = mod(iPlot - 1, obj.nCol) + 1;
                    iMark = mod(iPlot - 1, obj.nMarks) + 1;
                    pData = [Res(iN).steps.th];
                    ind = find(isnan(pData(1,:)),1,"first");
                    if isempty(ind); ind = steps(end) + 1; end
                    scatter(ax,max_steps,pData(ind),...
                        obj.marksSize * 1.5,obj.colorsFull(iCol,:),...
                        obj.marksOrder{iMark},...
                        LineWidth=obj.marksLW);
                end
            end
            %% Annotate
            xlabel(ax,'Iterations',...
                Interpreter='latex');
            ylabel(ax,'Parameter $\theta$',...
                Interpreter='latex');
            xlim(ax,[0 max_steps]);
            ylim(ax,[-pi/2 pi/2]);
            %% Finish plot
            Args2 = namedargs2cell(Args2);
            obj.postPlot(fg,ax,Args2{:});
        end
        function [fg,ax] = Energy(obj,V_range,Args)
            arguments
                obj
                V_range
                Args.Ebar
                Args.eps
                Args.rhoz
                Args.P
                Args.P_Tresh    =  1E-2
                Args.ymin       = -3
                Args.ymax       =  3
                Args.zooms      (1,:)   Plot.rect
                Args.zoomPos    (1,:)   Plot.rect
                Args.zoomTicks  = false
                Args.File
            end
            % ENERGY Plot the energy spectra $\theta$

            %% Initialize basic figures
            [fg,ax]=obj.prePlot;
            nLegends = 0;
            stateSuffix={'_+','_-'};
            %% Plot Data
            if isfield(Args,'P')
                %% Plot Energy spectra
                N = 2;
                k_max = (size(Args.P,3) - 1) / 2;
                for iN = 1:N
                    for k = -k_max:k_max
                        py = Args.eps(:,iN) + k * obj.calcObj.w;
                        % Skip if no energy is in window
                        ind = find(py > Args.ymin & py < Args.ymax);
                        if isempty(ind); continue; end
                        py = py(ind);
                        px = V_range(ind);
                        pmeta = Args.P(ind,iN,k_max+1+k);
                        % Skip if all states are barely visible
                        if max(pmeta) < Args.P_Tresh; continue; end
                        % Plot the spectra
		                patch(ax,[px,NaN],[py;NaN],obj.colors2{iN}, ...
                            'LineWidth',obj.lineWidth,'AlphaDataMapping','none',...
                            'EdgeAlpha','interp','EdgeColor',obj.colors2{iN},...
                            'FaceVertexAlphaData',[pmeta;nan]);
                    end
                    nLegends = nLegends + 1;
                    lgds(nLegends) = plot(ax,nan,nan,'-',...
                        'LineWidth',obj.lineWidth,'Color',obj.colors2{iN});
                    lgds_text{nLegends} = sprintf('$P^{(k)}%s$',stateSuffix{iN});
                end
            end
            if isfield(Args,'Ebar')
                %% Plot Average energies
                N = 2;
                for iN=1:N
                    plot(ax,V_range,Args.Ebar(:,iN),'--',...
                        'LineWidth',obj.lineWidth,'Color',obj.colors2{iN});
                    nLegends = nLegends + 1;
                    lgds(nLegends) = plot(ax,nan,nan,'--',...
                        'LineWidth',obj.lineWidth,'Color',obj.colors2{iN});
                    lgds_text{nLegends} = sprintf('$\\bar{E}%s$',stateSuffix{iN});
                end
            end
            %% Plot the zoomed-in section
            if isfield(Args,'zooms')
                if isfield(Args,'zoomPos')
                    if length(Args.zooms) ~= length(Args.zoomPos)
                        error('Size missmatch');
                    end
                else
                    if length(Args.zooms) > 3
                        error('Too many zoomed plots');
                    end
                    ax_pos = ax.Position;
                    ax_x0 = ax_pos(1);
                    ax_y0 = ax_pos(2);
                    ax_w = ax_pos(3);
                    ax_h = ax_pos(4);
                    y0 = ax_y0 + ax_h * 3/4 - ax_h / 20;
                    for iN = 1:length(Args.zooms)
                        x0 = ax_x0 + ax_w/20 + (iN - 1) * (ax_w/4+ax_w/10);
                        Args.zoomPos(iN) = Plot.rect(x0,y0,...
                            width=ax_w/4,height=ax_h/4);
                    end
                end
                for iN = 1:length(Args.zooms)
                    obj.zoomin(ax,Args.zooms(iN),Args.zoomPos(iN),...
                        Corners={'bot left','bot right'});
                end
            end
            %% Annotate
            switch length(lgds)
                case 4
                    lgds_cols = 2;
                otherwise
                    lgds_cols = 1;
            end
            legend(ax,lgds,lgds_text,...
                Location='south west',...
                Interpreter='latex',...
                NumColumns=lgds_cols);
            xlabel(ax,'Driving strength $V$',...
                Interpreter='latex');
            ylabel(ax,'Energy $E$',...
                Interpreter='latex');
            ylim(ax,[Args.ymin Args.ymax]);
            xlim(ax,V_range([1 end]));
            title(ax,'Two-level system without RWA');
            %% Finish plot
            subArgs = cell(1,0);
            if isfield(Args,'File')
                subArgs = [subArgs(:)',{'File',Args.File}];
            end
            obj.postPlot(fg,ax,subArgs{:});
        end
        function [fg,ax] = SteadyState(obj,V_range,Args)
            arguments
                obj
                V_range
                Args.rhoz
                Args.overlapM
                Args.File
            end
            % STEADYSTATE Plot the steady state distribution and overlap
            % measure

            %% Initialize basic figures
            [fg,ax]=obj.prePlot;
            %% Plot Data
            if isfield(Args,'rhoz')
                %% Plot the steady state
                if isfield(Args,'overlapM')
                    % Plot it with yyplot
                    yyaxis(ax,'left');
                end
                plot(ax,V_range,Args.rhoz,'-',...
                    LineWidth=obj.lineWidth,DisplayName='$\rho_z$');
            end
            if isfield(Args,'overlapM')
                %% Plot the steady state
                if isfield(Args,'rhoz')
                    % Plot it with yyplot
                    yyaxis(ax,'right');
                end
                plot(ax,V_range,Args.overlapM,'--',...
                    LineWidth=obj.lineWidth,DisplayName='Overlap Measure');
            end
            %% Annotate
%             legend(ax,...
%                 Location='south west',...
%                 Interpreter='latex');
            xlabel(ax,'Driving strength $V$',...
                Interpreter='latex');
            if isfield(Args,'rhoz')
                if isfield(Args,'overlapM')
                    yyaxis(ax,'left');
                end
                ylabel(ax,'SS density matrix $\rho_z$',...
                    Interpreter='latex');
                ylim(ax,[-1 1]);
            end
            if isfield(Args,'overlapM')
                if isfield(Args,'rhoz')
                    yyaxis(ax,'right');
                end
                ylabel(ax,'Overlap measure $\mathcal{S}$',...
                    Interpreter='latex');
                ylim(ax,[0 1]);
            end
            xlim(ax,V_range([1 end]));
            %% Finish plot
            subArgs = cell(1,0);
            if isfield(Args,'File')
                subArgs = [subArgs(:)',{'File',Args.File}];
            end
            obj.postPlot(fg,ax,subArgs{:});
        end
    end
end