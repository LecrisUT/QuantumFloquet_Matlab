classdef HubardDimer < handle
    methods
        function obj=HubardDimer()
        end
        function Res=Overlap(obj,U_range,t_range,Args)
            arguments
                obj
                U_range
                t_range
                Args.res_Exact
                Args.res_Orb
                Args.loadFile
                Args.saveFile
                Args.xi         = 1E-4
                Args.Print      = true
                Args.S_cutoff   = 1E-8
            end
            %% Calculate the non-SCF solutions
            nU = length(U_range);
            nt = length(t_range);
            if isfield(Args,'res_Exact')
                res_Exact = Args.res_Exact;
            end
            if isfield(Args,'loadFile')
                iU = load(Args.loadFile,'iU');
                Res = load(Args.loadFile,'Res');
                res_Exact = load(Args.loadFile,'res_Exact');
                res_Orb = load(Args.loadFile,'res_NSCF');
                NSol = size(res_Orb,2);
            elseif isfield(Args,'res_Orb')
                res_Orb = Args.res_Orb;
                NSol = size(res_Orb,2);
                Res(nU,NSol) = struct('S',[],'S_th',[],'t_max_th',[],'S_max',[],'S_min',[]);
                iU = 1;
            else
                error('No res_Orb deduced');
            end
            if isempty(whos('res_Exact'))
                % TODO: Perform Exact calculation if not provided
                error('Not implemented');
            end
            if isfield(Args,'saveFile')
                save(Args.saveFile,'Res','iU','res_Exact','res_Orb','obj','-v7.3');
            end
            %% Initialize missing variables
            k_max2 = length(res_Orb(1).Psi) / 2;
            k_max = (k_max2 - 1) / 2;
            %% Calculate the overlap
            for iU = iU:nU
                tic;
                tPsi_ex = reshape([res_Exact(iU,:).Psi],[],3);
                tPsit_ex = permute(reshape(tPsi_ex,3,k_max2,3),[1 3 2]);
                teps_ex = [res_Exact(iU,:).eps];
                for iN=1:NSol
                    if isempty(res_Orb(iU,iN).Psi)
                        % Populate nan for missing data
                        Res(iU,iN).S = nan(3,1);
                        Res(iU,iN).S_th = [];
                        Res(iU,iN).S_min = [];
                        Res(iU,iN).S_max = [];
                        Res(iU,iN).t_max_th = nan(4,1);
                    else
                        Res(iU,iN).S_th = zeros(4,nt);
                        Res(iU,iN).S_max = zeros(4,nt);
                        Res(iU,iN).S_min = zeros(4,nt);
                        %% Calculate the Floquet space overlap
                        tPsi = OrbToSlat(res_Orb(iU,iN).Psi);
                        tPsi0 = sum(tPsi,3);
                        S = zeros(3,1);
                        for jN=1:3
                            for k = -k_max:k_max
                                S(jN) = max(S(jN),...
                                    abs(circshift(tPsi_ex(:,jN),3*k,1)' * tPsi(:)));
                            end
                        end
                        Res(iU,iN).S = S;
                        %% Calculate the time-dependent overlaps
                        for it = 1:nt
                            t = t_range(it);
                            if it == 1
                                t_prev = 0;
                                prev_S_avg = zeros(4,1);
                            else
                                t_prev = t_range(it-1);
                                prev_S_avg = Res(iU,iN).S_th(:,it-1);
                            end
                            t_subrange=fliplr(t:-max_delt:t_prev);
                            tS_th = nan(4,1);
                            prev_S_max = nan(4,1);
                            prev_S_min = nan(4,1);
                            prev_t = t_prev;
                            for jt = 1:length(t_subrange)
                                tt = t_subrange(jt);
                                delt = tt - prev_t;
                                tS_th(1:3) = abs(PropFloquet(tt,w,tPsit_ex,'eps',teps_ex)' * ...
                                    PropFloquet(tt,w,tPsi));
                                tS_th(4) = abs(PropFloquet(tt,w,tPsit_ex,'eps',teps_ex,'Psi0',tPsi0)' * ...
                                    PropFloquet(tt,w,tPsi));
                                tS_min = min(prev_S_min,tS_th);
                                tS_max = max(prev_S_max,tS_th);
                                % Simple numerical integration
                                if prev_t < max_delt
                                    tS_avg = tS_th;
                                else
                                    tS_avg = (prev_S_avg * prev_t + tS_th * delt) / tt;
                                end
                                % If the overlaps do not change drasticly
                                if prod(tS_min == prev_S_min) && ...
                                        prod(tS_max == prev_S_max) && ...
                                        prod(abs(tS_avg - prev_S_avg) < Args.S_cutoff)
                                    tflag = true;
                                    tcount = tcount + 1;
                                    if tcount > 1
            %                             if tcount > 100000
                                        break
                                    end
                                else
                                    if tflag
                                        treset = treset+1;
            %                                 fprintf('Resetted [%d] %f\n',treset,log10(t));
                                    end
                                    tflag = false;
                                    tcount = 0;
                                end
                                prev_S_min = tS_min;
                                prev_S_max = tS_max;
                                prev_S_avg = tS_avg;
                                prev_t = tt;
                            end
                            Res(iU,iN).S_th(:,it) = tS_avg;
                            Res(iU,iN).S_min(:,it) = tS_min;
                            Res(iU,iN).S_max(:,it) = tS_max;
%                             fprintf('Done [%d,%d/%d,2] %E/%E\n',iU,iN,nU,t,t_range(end));
                        end
                    end
                    %% Find when overlap deviates
                    t_max = zeros(4,1);
                    tS = Res(iU,iN).S_th(:,1);
                    tS(4) = 1;
                    for ind = 1:4
                        tt_max = find(Res(iU,iN).S_min(ind,:) < tS(ind)-xi,1,'first');
                        if isempty(tt_max); tt_max = t_range(end); end
                        t_max(ind) = tt_max;
                    end
                    Res(iU,iN).t_max_th = t_max;
                end
                if isfield(Args,'saveFile')
                    save(Args.saveFile,'Res','iU','res_Exact','res_Orb','obj','-v7.3');
                end
                if Args.Print
                    fprintf('Done [%d/%d] (%fs)\n',iU,nU,toc);
                end
            end
        end
    end
end