clc;clear all;
w0=1; w=1.5;
del=w-w0;
sigma = [1 0 0];
system = QuanFloq.Models.TwoLevel(w0=w0,w=w,k_max=10,v=0,xi=1E-12,...
    RWA=true,rand_v=false);
bath_0K = QuanFloq.Models.BosonBath1(1,beta=inf);
bath_10B = QuanFloq.Models.BosonBath1(1,beta=10);
system_bath_0K = QuanFloq.Models.TwoLevel_Boson(sigma,system=system,bath=bath_0K);
system_bath_10B = QuanFloq.Models.TwoLevel_Boson(sigma,system=system,bath=bath_10B);
bath2_0K = QuanFloq.Models.BosonBath2(1,beta=inf);
bath2_10B = QuanFloq.Models.BosonBath2(1,beta=10);
system_bath2_0K = QuanFloq.Models.TwoLevel_Boson(sigma,system=system,bath=bath2_0K);
system_bath2_10B = QuanFloq.Models.TwoLevel_Boson(sigma,system=system,bath=bath2_10B);
fine_tune = [];
% fine_tune(1) = struct(Vc=sqrt(w^2-del^2)+1E-15, delV=4E-4, nV=1001,...
%     V_range=[]);
% fine_tune(2) = struct(Vc=sqrt(4*w^2-del^2)+1E-15, delV=4E-4, nV=1001,...
%     V_range=[]);
nV=501;
V_range=linspace(0,5,nV);
for iF = 1:length(fine_tune)
    ft = fine_tune(iF);
    fine_tune(iF).V_range = linspace(-ft.delV,ft.delV,ft.nV) + ft.Vc;
    ft = fine_tune(iF);
    V_range=unique(sort([V_range ft.V_range]));
end
nV = length(V_range);
%% Save calculation metadata
Details = QuanFloq.Details(script='Calc_TwoLevel3.m',...
    model='Two-level system',...
    calculation='Steady state calculation',...
    objects=[system bath_0K system_bath_0K bath_10B system_bath_10B ...
    bath2_0K system_bath2_0K bath2_10B system_bath2_10B],...
    variables=struct(V_range=V_range,fine_tune=fine_tune));
%% Calculate adiabatic eigenstates
Psi0 = [0 1;1 0];
eps0 = [-w0/2;w0/2];
iter=QuanFloq.GenericCalcIterator(system,data=V_range,updatefcn=@AdiabaticV);
iter.reset;
[tPsi,teps,tEbar] = system.eigs(iterator=iter,...
    Psi_prev=Psi0,eps_prev=eps0);
Res_system = PackRes(system,tPsi,teps,Ebar=tEbar);
%% Calculate steady-states
Res_0K(nV) = struct(rho=[],rhoz=[],OverlapMeasure=[],OverlapMeasure2=[]);
Res_10B(nV) = struct(rho=[],rhoz=[],OverlapMeasure=[]);
Res2_0K(nV) = struct(rho=[],rhoz=[],OverlapMeasure=[]);
Res2_10B(nV) = struct(rho=[],rhoz=[],OverlapMeasure=[]);
for iV = 1:nV
    system.V = V_range(iV);
    Psi = reshape([Res_system(iV,:).Psi],[],2);
    eps = [Res_system(iV,:).eps];
    Ebar = [Res_system(iV,:).Ebar];
    Overlap(1) = system.SpectraOverlap(Psi(:,1),Psi(:,2),eps1=eps(1),eps2=eps(2),sk_max=10);
    Overlap(2) = system.SpectraOverlap(Psi(:,2),Psi(:,1),eps1=eps(2),eps2=eps(1),sk_max=10);
    Overlap_0K(1) = system_bath_0K.SpectraOverlap(Psi(:,1),Psi(:,2),eps1=eps(1),eps2=eps(2),sk_max=10);
    Overlap_0K(2) = system_bath_0K.SpectraOverlap(Psi(:,2),Psi(:,1),eps1=eps(2),eps2=eps(1),sk_max=10);
    Overlap_10B(1) = system_bath_10B.SpectraOverlap(Psi(:,1),Psi(:,2),eps1=eps(1),eps2=eps(2),sk_max=10);
    Overlap_10B(2) = system_bath_10B.SpectraOverlap(Psi(:,2),Psi(:,1),eps1=eps(2),eps2=eps(1),sk_max=10);
    if Ebar(1) > Ebar(2)
        Overlap_0K(3) = Overlap_0K(2);
        Overlap_10B(3) = Overlap_10B(2);
        Overlap(3) = Overlap(2);
    else
        Overlap_0K(3) = Overlap_0K(1);
        Overlap_10B(3) = Overlap_10B(1);
        Overlap(3) = Overlap(1);
    end
    Overlap2_0K(1) = system_bath2_0K.SpectraOverlap(Psi1=Psi(:,1),Psi2=Psi(:,2),eps1=eps(1),eps2=eps(2),sk_max=10);
    Overlap2_0K(2) = system_bath2_0K.SpectraOverlap(Psi1=Psi(:,2),Psi2=Psi(:,1),eps1=eps(2),eps2=eps(1),sk_max=10);
    Overlap2_10B(1) = system_bath2_10B.SpectraOverlap(Psi1=Psi(:,1),Psi2=Psi(:,2),eps1=eps(1),eps2=eps(2),sk_max=10);
    Overlap2_10B(2) = system_bath2_10B.SpectraOverlap(Psi1=Psi(:,2),Psi2=Psi(:,1),eps1=eps(2),eps2=eps(1),sk_max=10);
    if Ebar(1) > Ebar(2)
        Overlap2_0K(3) = Overlap2_0K(2);
        Overlap2_10B(3) = Overlap2_10B(2);
    else
        Overlap2_0K(3) = Overlap2_0K(1);
        Overlap2_10B(3) = Overlap2_10B(1);
    end
    Gamma_0K = system_bath_0K.Gamma(Psi=Psi,eps=eps,type='Lindblad').Gamma;
    Gamma_10B = system_bath_10B.Gamma(Psi=Psi,eps=eps,type='Lindblad').Gamma;
    L_0K = system_bath_0K.Lindblad(Gamma=Gamma_0K);
    L_10B = system_bath_10B.Lindblad(Gamma=Gamma_10B);
    trho_0K = system_bath_0K.SteadyState(type='Lindblad',MEq=L_0K);
    trho_10B = system_bath_10B.SteadyState(type='Lindblad',MEq=L_10B);
    if isempty(trho_0K) || isempty(trho_10B)
        error('No steady states');
    end
    if size(trho_0K,2) > 1 || size(trho_10B,2) > 1
        warning('Multiple steady states');
    end
    if ~isempty(find(trho_0K<0,1)) || ~isempty(find(trho_10B<0,1))
        error('Negative population');
    end
    Gamma2_0K = system_bath2_0K.Gamma(Psi=Psi,eps=eps,type='Lindblad').Gamma;
    Gamma2_10B = system_bath2_10B.Gamma(Psi=Psi,eps=eps,type='Lindblad').Gamma;
    L2_0K = system_bath_0K.Lindblad(Gamma=Gamma2_0K);
    L2_10B = system_bath_10B.Lindblad(Gamma=Gamma2_10B);
    trho2_0K = system_bath_0K.SteadyState(type='Lindblad',MEq=L2_0K);
    trho2_10B = system_bath_10B.SteadyState(type='Lindblad',MEq=L2_10B);
    if isempty(trho2_0K) || isempty(trho2_10B)
        error('No steady states');
    end
    if size(trho2_0K,2) > 1 || size(trho2_10B,2) > 1
        warning('Multiple steady states');
    end
    if ~isempty(find(trho2_0K<0,1)) || ~isempty(find(trho2_10B<0,1))
        error('Negative population');
    end
    Res_0K(iV).rho = trho_0K;
    Res_0K(iV).rhoz = trho_0K(2) - trho_0K(1);
    Res_0K(iV).OverlapMeasure = Overlap_0K;
    Res_0K(iV).OverlapMeasure2 = Overlap(3);
    Res_10B(iV).rho = trho_10B;
    Res_10B(iV).rhoz = trho_10B(2) - trho_10B(1);
    Res_10B(iV).OverlapMeasure = Overlap_10B;
    Res2_0K(iV).rho = trho2_0K;
    Res2_0K(iV).rhoz = trho2_0K(2) - trho2_0K(1);
    Res2_0K(iV).OverlapMeasure = Overlap2_0K;
    Res2_10B(iV).rho = trho2_10B;
    Res2_10B(iV).rhoz = trho2_10B(2) - trho2_10B(1);
    Res2_10B(iV).OverlapMeasure = Overlap2_10B;
    for iN = 1:2
        Res_system(iV,iN).OverlapMeasure = Overlap(iN);
    end
    fprintf('Done %d/%d\n',iV,nV);
end
%% PostProcess
%% Save results
save('Calc_TwoLevel3.mat',...
    'Res_system','Res_0K','Res_10B',...
    'Res2_0K','Res2_10B',...
    'Details',...
    '-v7.3');
%% Helper functions
function AdiabaticV(obj)
    obj.object.V = obj.data(obj.ind);
end
function Res = PackRes(obj,Psi,eps,Args)
    arguments
        obj
        Psi
        eps
        Args.Ebar
    end
    M = size(Psi,2);
    L = size(Psi,3);
    Res(L,M) = struct(Psi=[],eps=[],P=[]);
    for iL = 1:L
        tPsi = obj.Psi_Fourier(Psi(:,:,iL));
        for iM = 1:M
            %% Pack calculated results
            Res(iL,iM).Psi = Psi(:,iM,iL);
            Res(iL,iM).eps = eps(iM,iL);
            if isfield(Args,'Ebar')
                Res(iL,iM).Ebar = Args.Ebar(iM,iL);
            end
            %% Other PostProcess
            Res(iL,iM).P = zeros(obj.k_max2,1);
            for ik = 1:obj.k_max2
                Res(iL,iM).P(ik) = tPsi(:,iM,ik)' * tPsi(:,iM,ik);
            end
        end
    end
end