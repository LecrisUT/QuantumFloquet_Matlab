clc; clear all;
w=8.3; N = 20;
xi=1E-2; tol=1E-12;
obj = QuanFloq.Models.BoxParticle(N,w=w,k_max=200,xi=1E-12);
fine_tune(1) = struct(Vc=48.97, delV=1E-1, nV=41,...
    V_range=[]);
fine_tune(2) = struct(Vc=34.08, delV=1E-1, nV=41,...
    V_range=[]);
fine_tune(3) = struct(Vc=45.8, delV=2E-1, nV=41,...
    V_range=[]);
nV=501;
V_range=linspace(0,50,nV);
for iF = 1:length(fine_tune)
    ft = fine_tune(iF);
    fine_tune(iF).V_range = linspace(-ft.delV,ft.delV,ft.nV) + ft.Vc;
    ft = fine_tune(iF);
    V_range=unique(sort([V_range ft.V_range]));
end
nV = length(V_range);
%% Save calculation metadata
Details = QuanFloq.Details(script='Calc_BoxParticle.m',...
    model='Particle in a box',...
    calculation='Adiabatic Floquet eigenstate calculation',...
    objects=obj,...
    variables=struct(V_range=V_range,fine_tune=fine_tune));
%% Calculate eigenstates
% Adiabatically continue from static eigenstates
Psi0 = eye(N);
eps0 = (1:N).^2;
iter=QuanFloq.GenericCalcIterator(obj,data=V_range,updatefcn=@AdiabaticV);
iter.reset;
[tPsi,teps,tEbar] = obj.eigs(iterator=iter,...
    Psi_prev=Psi0,eps_prev=eps0);
Res_AE = PackRes(obj,tPsi,teps,Ebar=tEbar);
iter.reset;
[tPsi,teps] = obj.eigs(iterator=iter,...
    Psi_prev=Psi0,eps_prev=eps0);
Res_QE = PackRes(obj,tPsi,teps);
%% PostProcess
%% Save results
save('Calc_BoxParticle.mat',...
    'Res_QE','Res_AE',...
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