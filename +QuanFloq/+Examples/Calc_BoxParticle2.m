clc; clear all;
N=13;
w=8.3;
xi=1E-6;
tol=1E-12;
k_max = 15;
obj = QuanFloq.Models.BoxParticle(N,w=w,k_max=k_max,xi=xi);
obj_ex = QuanFloq.Models.BoxParticle(20,w=w,k_max=200,xi=tol);
V = 30;
obj.V = V;
obj_ex.V = V;
NRand = 40;
% Psi0=1-2*rand(N,NRand);
Psi0=zeros(N,NRand);
Psi0(1:7,:)=1-2*rand(7,NRand);
% Psi0(1:10,:)=1-2*rand(10,NRand);
Psi0=Psi0./vecnorm(Psi0);

%% Save calculation metadata
Details = QuanFloq.Details(script='Calc_BoxParticle2.m',...
    model='Particle in a box',...
    calculation='Variational Floquet eigenstate convergence',...
    objects=[obj,obj_ex],...
    variables=struct(V=V,Psi0=Psi0));
%% Calculate exact eigenstates
[Res_Ex.Psi,Res_Ex.eps,Res_Ex.Ebar] = obj_ex.eigs;
%% Calculate eigenstates
%     Display='iter',...
warning('off', 'MATLAB:nearlySingularMatrix');
warning('off', 'MATLAB:SingularMatrix');
%     Display='iter',...
opts=optimoptions('fmincon',Algorithm='interior-point',...
    FunctionTolerance=xi,OptimalityTolerance=xi,...
    StepTolerance=tol,ConstraintTolerance=tol,...
    MaxIterations=1E2,MaxFunctionEvaluations=1E5);
fprintf('(Unpert) Doing AE variational principle:\n');tic;
Res_AE = obj.variational('AE',Psi0=Psi0,...
    TrackPsi=true,filter_nconv=false,tol=tol,...
    opts=opts);
fprintf('(Unpert) Done AE method: %f (s)\n',toc);
% opts=optimoptions(opts,...
%     FunctionTolerance=tol,OptimalityTolerance=tol,...
%     StepTolerance=tol,ConstraintTolerance=tol);
% fprintf('(Unpert) Doing QE variance variational principle:\n');tic;
% Res_QE = obj.variational('varQE',Psi0=Psi0,...
%     TrackPsi=true,filter_nconv=false,tol=tol,...
%     opts=opts);
% fprintf('(Unpert) Done QE method: %f (s)\n',toc);
warning('on', 'MATLAB:nearlySingularMatrix');
warning('on', 'MATLAB:SingularMatrix');
%% PostProcess
Processer = Process.BoxParticle(obj,exactObj=obj_ex);
Res_AE = Processer.Variational(Res_AE,...
    Norm=true,Ex_Psi=Res_Ex.Psi,...
    Ex_varEps=true);
% Res_QE = Processer.Variational(Res_QE,...
%     Norm=true,Ex_Psi=Res_Ex.Psi,...
%     Ex_varEps=true);
%% Save results
save('Calc_BoxParticle2.mat',...
    'Res_AE',...
    'Details',...
    '-v7.3');
% save('Calc_BoxParticle2.mat',...
%     'Res_QE','Res_AE',...
%     'Details',...
%     '-v7.3');
