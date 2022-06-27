clc; clear all;
w0=1; w=1.5;
del=w-w0;
xi=1E-2; tol=1E-12; v=1E-4;
Unpert = Calc.TwoLevel(w0=w0,w=w,k_max=100,v=0,xi=tol);
Pert = Calc.TwoLevel(w0=w0,w=w,k_max=100,v=v,xi=xi);
Pert2 = Calc.TwoLevel(w0=w0,w=w,k_max=100,v=v,xi=tol);
V = sqrt(w^2-del^2);
Unpert.V = V;
Pert.V = V;
Pert2.V = V;
th_range=linspace(-pi/2,pi/2,21);
nth = length(th_range);
%% Save calculation metadata
Details = Calc.Details(script='Calc_TwoLevel2.m',...
    model='Two-level system',...
    calculation='Variational Floquet eigenstate convergence',...
    objects=[Pert,Unpert,Pert2],...
    variables=struct(V=V,th_range=th_range));
%% Calculate exact eigenstates
[Res_Ex.eps,Res_Ex.Psi,Res_Ex.Ebar] = Unpert.eigs;
%% Calculate eigenstates
fprintf('(Unpert) Doing AE variational principle:\n');tic;
Res_Unpert_AE = Unpert.variational('AE',th_range=th_range,...
    TrackPsi=true,filter_nconv=false,tol=tol);
fprintf('(Unpert) Done AE method: %f (s)\n',toc);
fprintf('(Unpert) Doing QE variance variational principle:\n');tic;
Res_Unpert_QE = Unpert.variational('varQE',th_range=th_range,...
    TrackPsi=true,filter_nconv=false,tol=tol);
fprintf('(Unpert) Done QE method: %f (s)\n',toc);
fprintf('(Pert) Doing AE variational principle:\n');tic;
Res_Pert_AE = Pert.variational('AE',th_range=th_range,...
    TrackPsi=true,filter_nconv=false,tol=tol);
fprintf('(Pert) Done AE method: %f (s)\n',toc);
fprintf('(Pert) Doing QE variance variational principle:\n');tic;
Res_Pert_QE = Pert.variational('varQE',th_range=th_range,...
    TrackPsi=true,filter_nconv=false,tol=tol);
fprintf('(Pert) Done QE method: %f (s)\n',toc);
fprintf('(Pert) Doing AE variational principle:\n');tic;
Res_Pert2_AE = Pert2.variational('AE',th_range=th_range,...
    TrackPsi=true,filter_nconv=false,tol=tol);
fprintf('(Pert) Done AE method: %f (s)\n',toc);
%% PostProcess
Processer = Process.TwoLevel(Unpert);
Res_Unpert_AE = Processer.Variational(Res_Unpert_AE,...
    th=true,Norm=true,Ex_Psi=Res_Ex.Psi);
Res_Unpert_QE = Processer.Variational(Res_Unpert_QE,...
    th=true,Norm=true,Ex_Psi=Res_Ex.Psi);
Processer.calcObj = Pert;
Res_Pert_AE = Processer.Variational(Res_Pert_AE,...
    th=true,Norm=true,Ex_Psi=Res_Ex.Psi);
Res_Pert_QE = Processer.Variational(Res_Pert_QE,...
    th=true,Norm=true,Ex_Psi=Res_Ex.Psi);
Processer.calcObj = Pert2;
Res_Pert2_AE = Processer.Variational(Res_Pert2_AE,...
    th=true,Norm=true,Ex_Psi=Res_Ex.Psi);
%% Save results
save('Calc_TwoLevel2.mat',...
    'Res_Unpert_QE','Res_Unpert_AE',...
    'Res_Pert_QE','Res_Pert_AE','Res_Pert2_AE',...
    'Details',...
    '-v7.3');

