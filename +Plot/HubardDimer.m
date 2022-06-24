classdef HubardDimer < Plot.basePlot
    properties
        colors3     = {'red','green','blue'}
        lineOrder   = {'-.','--'}
    end
    properties (Dependent)
        nLines
    end
    methods
        function val = get.nLines(obj)
            val = length(obj.lineOrder);
        end
        function obj=HubardDimer(calcObj)
            obj@Plot.basePlot(calcObj);
        end
        function PlotAll(obj,Args)
            arguments
                obj
                Args.Exact
                Args.FHF
                Args.NSCF
            end
            % PLOTALL Plot all default

            
        end
        function Theta(obj,Args)
            arguments
                obj
                Args.Exact
                Args.FHF
                Args.NSCF
                Args.File
            end
            % THETA Plot the parameter $\theta$

            %% Initialize basic figures
            [fg,ax]=obj.prePlot;
            nLegends = 0;
            xmin=nan; xmax=nan;
            %% Plot Data
            if isfield(Args,'Exact')
                %% Plot Exact solutions data
                U_range = Args.Exact.U_range;
                xmin=min(xmin,U_range(1));
                xmax=min(xmax,U_range(end));
                nU = length(U_range);
                N = 3;
                NBranch = 2;
                pData = permute(reshape([Args.Exact.Res(:,:).th],...
                    NBranch,nU,N),[2 3 1]);
                pMeta = permute(reshape([Args.Exact.Res(:,:).S],...
                    NBranch,nU,N),[2 3 1]);
                [U_range,pData,pMeta]=obj.breakLines(U_range,pData,'meta',pMeta);
                for iN = 1:N
                    for iB = 1:NBranch
                        patch(ax,[U_range nan],[pData(:,iN,iB);nan],obj.colors3{iN},...
                            'LineWidth',obj.lineWidth,'AlphaDataMapping','none',...
                            'EdgeAlpha','interp','EdgeColor',obj.colors3{iN},...
                            'FaceVertexAlphaData',[pMeta(:,iN,iB);nan]);
                    end
                end
                nLegends = nLegends + 1;
                lgds(nLegends) = plot(ax,nan,nan,'-',...
                    'LineWidth',obj.lineWidth,'Color',obj.colors3{1});
                lgds_text{nLegends} = '$\Phi$ Exact';
            end
            if isfield(Args,'FHF')
                %% Plot Floquet Hartree-Fock data
                U_range = Args.FHF.U_range;
                xmin=min(xmin,U_range(1));
                xmax=min(xmax,U_range(end));
                nU = length(U_range);
                N = Args.FHF.NSol;
                pData = permute(reshape([Args.FHF.Res(:,:).th],...
                    1,nU,N),[2 3 1]);
                for iN=1:N
                    iCol = mod(iN - 1, obj.nCol) + 1;
                    iMark = mod(iN - 1, obj.nMarks) + 1;
                    plot(ax,U_range,pData(:,iN),...
                        'LineWidth',obj.marksLW,'MarkerSize',obj.marksSize,...
                        'LineStyle','none','Color',obj.colorsFull(iCol,:),...
                        'Marker',obj.marksOrder{iMark});
                end
                nLegends = nLegends + 1;
                lgds(nLegends) = plot(ax,nan,nan,obj.marksOrder{1},...
                    'LineWidth',obj.marksLW,'Color',obj.colorsFull(iCol,:));
                lgds_text{nLegends} = '$\phi$ Floquet HF';
            end
            if isfield(Args,'NSCF')
                %% Plot non-SCF data
                U_range = Args.NSCF.U_range;
                xmin=min(xmin,U_range(1));
                xmax=min(xmax,U_range(end));
                nU = length(U_range);
                N = 2;
                pData = permute(reshape([Args.NSCF.Res(:,:).th],...
                    1,nU,N),[2 3 1]);
                [U_range,pData]=obj.breakLines(U_range,pData);
                pColor={obj.colors3{1} obj.colors3{3}};
                for iN=1:N
                    plot(ax,U_range,pData(:,iN),obj.lineOrder{iN},...
                        'LineWidth',obj.lineWidth,'Color',pColor{iN});
                end
                nLegends = nLegends + 1;
                lgds(nLegends) = plot(ax,nan,nan,obj.lineOrder{1},...
                    'LineWidth',obj.lineWidth,'Color',pColor{1});
                lgds_text{nLegends} = 'NSCF';
            end
            %% Annotate
            legend(ax,lgds,lgds_text,...
                'Location','north west',...
                'Interpreter','latex');
            xlabel(ax,'Cuolomb interaction $U$',...
                'Interpreter','latex');
            ylabel(ax,'Parameter $\theta$',...
                'Interpreter','latex');
            ylim(ax,[-pi/2 pi/2]);
            xlim(ax,[xmin xmax]);
            title(ax,'Slater determinant closeness');
            %% Finish plot
            subArgs = cell(1,0);
            if isfield(Args,'File')
                subArgs = [subArgs(:)',{'File'},{Args.File}];
            end
            obj.postPlot(fg,ax,subArgs{:});
        end
        function [fg,ax]=AverageEnergy(obj,Args)
            arguments
                obj
                Args.Exact
                Args.FHF
                Args.File
            end
            % AVERAGEENERGY Plot the average energies

            %% Initialize basic figures
            [fg,ax]=obj.prePlot;
            nLegends = 0;
            xmin=nan; xmax=nan;
            %% Plot Data
            if isfield(Args,'Exact')
                %% Plot Exact solutions data
                U_range = Args.Exact.U_range;
                xmin=min(xmin,U_range(1));
                xmax=min(xmax,U_range(end));
                nU = length(U_range);
                N = 3;
                pData = permute(reshape([Args.Exact.Res(:,:).Ebar],...
                    1,nU,N),[2 3 1]);
                for iN = 1:N
                    plot(ax,U_range,pData(:,iN),'-',...
                        'LineWidth',obj.lineWidth,...
                        'Color',obj.colors3{iN});
                end
                nLegends = nLegends + 1;
                lgds(nLegends) = plot(ax,nan,nan,'-',...
                    'LineWidth',obj.lineWidth,'Color',obj.colors3{1});
                lgds_text{nLegends} = '$\Phi$ Exact';
            end
            if isfield(Args,'FHF')
                %% Plot Floquet Hartree-Fock data
                U_range = Args.FHF.U_range;
                xmin=min(xmin,U_range(1));
                xmax=min(xmax,U_range(end));
                nU = length(U_range);
                N = Args.FHF.NSol;
                pData = permute(reshape([Args.FHF.Res(:,:).Ebar],...
                    1,nU,N),[2 3 1]);
                for iN=1:N
                    iCol = mod(iN - 1, obj.nCol) + 1;
                    iMark = mod(iN - 1, obj.nMarks) + 1;
                    plot(ax,U_range,pData(:,iN),...
                        'LineWidth',obj.marksLW,'MarkerSize',obj.marksSize,...
                        'LineStyle','none','Color',obj.colorsFull(iCol,:),...
                        'Marker',obj.marksOrder{iMark});
                end
                nLegends = nLegends + 1;
                lgds(nLegends) = plot(ax,nan,nan,obj.marksOrder{1},...
                    'LineWidth',obj.marksLW,'Color',obj.colorsFull(1,:));
                lgds_text{nLegends} = '$\phi$ Floquet HF';
            end
            %% Annotate
            legend(ax,lgds,lgds_text,...
                'Location','north west',...
                'Interpreter','latex');
            xlabel(ax,'Cuolomb interaction $U$',...
                'Interpreter','latex');
            ylabel(ax,'Average Energy $\bar{E}$',...
                'Interpreter','latex');
%             ylim(ax,[-3 5]);
            xlim(ax,[xmin xmax]);
            title(ax,'Average energy');
            %% Finish plot
            subArgs = cell(1,0);
            if isfield(Args,'File')
                subArgs = [subArgs(:)',{'File'},{Args.File}];
            end
            obj.postPlot(fg,ax,subArgs{:});
        end
        function [fgs,axs]=OverlapEigen(obj,NEigen,Args)
            arguments
                obj
                NEigen      = 1:3
                Args.merged = false
                Args.NSCF
                Args.FHF
                Args.Files
            end
            % OVERLAPEIGEN Plot the overlap with exact eigenstates
            if Args.merged
                [fg,ax]=obj.prePlot;
                nPlots = 1;
            else
                nPlots = 0;
            end
            for iEigen = NEigen
                % Loop over each eigenstate
                %% Initialize basic figures
                if ~Args.merged
                    [fg,ax]=obj.prePlot;
                    nPlots = nPlots + 1;
                else
                    hold(ax,'on');
                end
                nLegends = 0;
                xmin=nan; xmax=nan;
                %% Plot Data
                if isfield(Args,'NSCF')
                    %% Plot non-SCF data
                    U_range = Args.NSCF.U_range;
                    xmin=min(xmin,U_range(1));
                    xmax=min(xmax,U_range(end));
                    nU = length(U_range);
                    N = 2;
                    NBranch = 3;
                    pData = permute(reshape([Args.NSCF.Res(:,:).S],...
                        NBranch,nU,N),[2 3 1]);
                    pData = pData(:,:,iEigen);
                    for iN=1:N
                        plot(ax,U_range,pData(:,iN),obj.lineOrder{iN},...
                            'LineWidth',obj.lineWidth,'Color',obj.colors3{iEigen});
                    end
                    nLegends = nLegends + 1;
                    lgds(nLegends) = plot(ax,nan,nan,obj.lineOrder{1},...
                        'LineWidth',obj.lineWidth,'Color',obj.colors3{iEigen});
                    lgds_text{nLegends} = 'NSCF';
                end
                if isfield(Args,'FHF')
                    %% Plot Floquet Hartree-Fock data
                    U_range = Args.FHF.U_range;
                    xmin=min(xmin,U_range(1));
                    xmax=min(xmax,U_range(end));
                    nU = length(U_range);
                    N = Args.FHF.NSol;
                    NBranch = 3;
                    pData = permute(reshape([Args.FHF.Res(:,:).S],...
                        NBranch,nU,N),[2 3 1]);
                    pData = pData(:,:,iEigen);
                    for iN=1:N
                        iMark = mod(iN - 1, obj.nMarks) + 1;
                        plot(ax,U_range,pData(:,iN),...
                            'LineWidth',obj.marksLW,'MarkerSize',obj.marksSize,...
                            'LineStyle','none','Color',obj.colors3{iEigen},...
                            'Marker',obj.marksOrder{iMark});
                    end
                    nLegends = nLegends + 1;
                    lgds(nLegends) = plot(ax,nan,nan,obj.marksOrder{1},...
                        'LineWidth',obj.marksLW,'Color',obj.colors3{iEigen});
                    lgds_text{nLegends} = '$\phi$ Floquet HF';
                end
                %% Annotate
                legend(ax,lgds,lgds_text,...
                    'Location','north west',...
                    'Interpreter','latex');
                xlabel(ax,'Cuolomb interaction $U$',...
                    'Interpreter','latex');
                ylabel(ax,'Overlap',...
                    'Interpreter','latex');
                ylim(ax,[0 1.1]);
                xlim(ax,[xmin xmax]);
                if Args.merged
                    ptitle = 'Overlap with exact eigenstates';
                else
                    ptitle = sprintf('Overlap with exact eigenstate $\\Psi_%d$',iEigen-1);
                end
                title(ax,ptitle,...
                    'Interpreter','latex');
                %% Finish plot
                % Make File field to save to
                subArgs = cell(1,0);
                if isfield(Args,'Files')
                    subArgs = [subArgs(:)',{'File',Args.Files{iEigen}}];
                end
                obj.postPlot(fg,ax,subArgs{:});
                fgs(nPlots) = fg;
                axs(nPlots) = ax;
            end
        end
        function [fgs,axs]=OverlapProp(obj,NEigen,Args)
            arguments
                obj
                NEigen          = 1:3
                Args.NSCF
                Args.FHF
                Args.Files
                Args.plot_S_th  = true
                Args.plot_S_min = true
            end
            % OVERLAPPROP Plot the propagated overlaps with the exact eigenstates

            fileInd = 0;
            %% Plot Data
            if isfield(Args,'NSCF')
                %% Prepare data
                U_range = Args.NSCF.U_range;
                nU = length(U_range);
                t_range = Args.NSCF.t_range;
                N = 2;
                for iN = 1:N
                    %% Get the current data
                    subArgs = cell(1,0);
                    if Args.plot_S_th
                        pData = cell(1,nU);
                        for iU = 1:nU
                            pData{iU} = Args.NSCF.Res(iU,iN).S_th;
                        end
                        subArgs = [subArgs(:)', {'S_th',pData}];
                    end
                    if Args.plot_S_min
                        pData1 = cell(1,nU);
                        pData2 = cell(1,nU);
                        for iU = 1:nU
                            pData1{iU} = Args.NSCF.Res(iU,iN).S_min;
                            pData2{iU} = Args.NSCF.Res(iU,iN).S_max;
                        end
                        subArgs = [subArgs(:)', {'S_min',pData1,'S_max',pData2}];
                    end
                    %% Plot non-SCF data
                    [fgs,axs]=obj.OverlapPropU(U_range,t_range,NEigen,subArgs{:});
                    %% Adjust plot

                    for iax = 1:length(axs)
                        ptitle = sprintf('NSCF overlap of linestyle [%s]\n',obj.lineOrder{iN});
                        axs(iax).Title.String = [ptitle axs(iax).Title.String];
                        %% Finish plot
                        % Make File field to save to
                        subArgs = cell(1,0);
                        if isfield(Args,'Files')
                            fileInd = fileInd + 1;
                            subArgs = [subArgs(:)',{'File',Args.Files{fileInd}}];
                        end
                        obj.postPlot(fgs(iax),axs(iax),subArgs{:});
                    end
                end
            end
            if isfield(Args,'FHF')
                %% Prepare data
                U_range = Args.FHF.U_range;
                nU = length(U_range);
                t_range = Args.FHF.t_range;
                N = Args.FHF.NSol;
                for iN = 1:N
                    %% Get the current data
                    subArgs = cell(1,0);
                    if Args.plot_S_th
                        pData = cell(1,nU);
                        for iU = 1:nU
                            pData{iU} = Args.FHF.Res(iU,iN).S_th;
                        end
                        subArgs = [subArgs(:)', {'S_th',pData}];
                    end
                    if Args.plot_S_min
                        pData1 = cell(1,nU);
                        pData2 = cell(1,nU);
                        for iU = 1:nU
                            pData1{iU} = Args.FHF.Res(iU,iN).S_min;
                            pData2{iU} = Args.FHF.Res(iU,iN).S_max;
                        end
                        subArgs = [subArgs(:)', {'S_min',pData1,'S_max',pData2}];
                    end
                    %% Plot Floquet HF data
                    [fgs,axs]=obj.OverlapPropU(U_range,t_range,NEigen,subArgs{:});
                    %% Adjust plot

                    for iax = 1:length(axs)
                        ptitle = sprintf('Floquet HF overlap of marker [%s]\n',obj.marksLatex{iN});
                        axs(iax).Title.String = [ptitle axs(iax).Title.String];
                        %% Finish plot
                        % Make File field to save to
                        subArgs = cell(1,0);
                        if isfield(Args,'Files')
                            fileInd = fileInd + 1;
                            subArgs = [subArgs(:)',{'File',Args.Files{fileInd}}];
                        end
                        obj.postPlot(fgs(iax),axs(iax),subArgs{:});
                    end
                end
            end
        end
        function [fgs,axs]=OverlapPropU(obj,U_range,t_range,NEigen,Args)
            arguments
                obj
                U_range
                t_range
                NEigen      = 1:4
                Args.S_th
                Args.S_min
                Args.S_max
                Args.Files
            end
            % OVERLAPPROPU Plot the propagated overlap with the eigenstates

            nPlots = 0;
            for iEigen = NEigen
                % Loop over each eigenstate
                nPlots = nPlots + 1;
                %% Initialize basic figures
                [fg,ax]=obj.prePlot;
                %% Prepare data
                nU = length(U_range);
                % Colorbar
                cbar=jet(nU);
                %% Plot Data
                colormap(ax,cbar);
                for iU=1:nU
                    if isfield(Args,'S_max')
                        % Plot contour of overlap range
                        if ~isempty(Args.S_max{iU})
                            patch(ax,[t_range fliplr(t_range)],...
                                [Args.S_min{iU}(iEigen,:) fliplr(Args.S_max{iU}(iEigen,:))],...
                                cbar(iU,:),FaceAlpha=0.4,EdgeColor=cbar(iU,:),EdgeAlpha=0.4);
                        end
                    end
                    if isfield(Args,'S_th')
                        % Plot Average overlap
                        if ~isempty(Args.S_th{iU})
                            semilogx(ax,t_range,Args.S_th{iU}(iEigen,:),...
                                Color=cbar(iU,:));
                        end
                    end
                end
                %% Annotate
                xlabel(ax,'Time $t$',...
                    'Interpreter','latex');
                ylabel(ax,'Overlap',...
                    'Interpreter','latex');
                if isfield(Args,'S_th')
                    ptitle = 'Time-averaged overlap';
                else
                    ptitle = 'Overlap boundaries';
                end
                if iEigen == 4
                    ptitle = sprintf('%s with Exact Slater determinant $\\Psi_{\\theta}$',ptitle);
                else
                    ptitle = sprintf('%s with Exact eigenstate $\\Psi_{%d}$',...
                        ptitle,iEigen-1);
                end
                title(ax,ptitle,...
                    'Interpreter','latex');
                caxis(ax,U_range([1 end]));
                ax.XScale="log";
                colorbar(ax);
                ylim(ax,[0 1.1]);
                %% Finish plot
                % Make File field to save to
                subArgs = cell(1,0);
                if isfield(Args,'Files')
                    subArgs = [subArgs(:)',{'File',Args.Files{iEigen}}];
                end
                obj.postPlot(fg,ax,subArgs{:});
                fgs(nPlots) = fg;
                axs(nPlots) = ax;
            end
        end
        function [fgs,axs]=t_max(obj,NEigen,Args)
            arguments
                obj
                NEigen      = 1:4
                Args.merged = false
                Args.NSCF
                Args.FHF
                Args.Files
            end
            % OVERLAPEIGEN Plot the overlap with exact eigenstates
            
            if Args.merged
                [fg,ax]=obj.prePlot;
                nPlots = 1;
            else
                nPlots = 0;
            end
            for iEigen = NEigen
                % Loop over each eigenstate
                %% Initialize basic figures
                if ~Args.merged
                    [fg,ax]=obj.prePlot;
                    nPlots = nPlots + 1;
                else
                    hold(ax,'on');
                end
                nLegends = 0;
                xmin=nan; xmax=nan;
                %% Plot Data
                if isfield(Args,'NSCF')
                    %% Plot non-SCF data
                    U_range = Args.NSCF.U_range;
                    xmin=min(xmin,U_range(1));
                    xmax=min(xmax,U_range(end));
                    nU = length(U_range);
                    N = 2;
                    NBranch = 4;
                    pData = permute(reshape([Args.NSCF.Res(:,:).t_max_th],...
                        NBranch,nU,N),[2 3 1]);
                    pData = pData(:,:,iEigen);
                    pColor = cell(1,N);
                    for iN=1:N
                        if iEigen == 4
                            pColor={obj.colors3{1} obj.colors3{3}};
                        else
                            pColor{iN} = obj.colors3{iEigen};
                        end
                        plot(ax,U_range,pData(:,iN),obj.lineOrder{iN},...
                            'LineWidth',obj.lineWidth,'Color',pColor{iN});
                    end
                    nLegends = nLegends + 1;
                    lgds(nLegends) = plot(ax,nan,nan,obj.lineOrder{1},...
                        'LineWidth',obj.lineWidth,'Color',pColor{1});
                    lgds_text{nLegends} = 'NSCF';
                end
                if isfield(Args,'FHF')
                    %% Plot Floquet Hartree-Fock data
                    U_range = Args.FHF.U_range;
                    xmin=min(xmin,U_range(1));
                    xmax=min(xmax,U_range(end));
                    nU = length(U_range);
                    N = Args.FHF.NSol;
                    NBranch = 4;
                    pData = permute(reshape([Args.FHF.Res(:,:).t_max_th],...
                        NBranch,nU,N),[2 3 1]);
                    pData = pData(:,:,iEigen);
                    pColor = cell(1,N);
                    for iN=1:N
                        if iEigen == 4
                            iCol = mod(iN - 1, obj.nCol) + 1;
                            pColor{iN} = obj.colorsFull(iCol,:);
                        else
                            pColor{iN} = obj.colors3{iEigen};
                        end
                        iMark = mod(iN - 1, obj.nMarks) + 1;
                        plot(ax,U_range,pData(:,iN),...
                            'LineWidth',obj.marksLW,'MarkerSize',obj.marksSize,...
                            'LineStyle','none','Color',pColor{iN},...
                            'Marker',obj.marksOrder{iMark});
                    end
                    nLegends = nLegends + 1;
                    lgds(nLegends) = plot(ax,nan,nan,obj.marksOrder{1},...
                        'LineWidth',obj.marksLW,'Color',pColor{1});
                    lgds_text{nLegends} = '$\phi$ Floquet HF';
                end
                %% Annotate
                legend(ax,lgds,lgds_text,...
                    'Location','north west',...
                    'Interpreter','latex');
                xlabel(ax,'Cuolomb interaction $U$',...
                    'Interpreter','latex');
                ylabel(ax,'Divergence point $T_max$',...
                    'Interpreter','latex');
%                 ylim(ax,[0 1.1]);
                xlim(ax,[xmin xmax]);
                ax.YScale="log";
                if Args.merged
                    ptitle = 'First divergence points from exact eigenstates';
                else
                    if iEigen == 4
                        ptitle='First divergence point from exact propagator $\Psi_\theta$';
                    else
                        ptitle=sprintf('First divergence point from eigenstate $\\Psi_%d$',iEigen-1);
                    end
                end
                title(ax,ptitle,...
                    'Interpreter','latex');
                %% Finish plot
                % Make File field to save to
                subArgs = cell(1,0);
                if isfield(Args,'Files')
                    subArgs = [subArgs(:)',{'File',Args.Files{iEigen}}];
                end
                obj.postPlot(fg,ax,subArgs{:});
                fgs(nPlots) = fg;
                axs(nPlots) = ax;
            end
        end
    end
end