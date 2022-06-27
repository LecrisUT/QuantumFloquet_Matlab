clc; clear all;
w0=1; w=1.5;
del=w-w0;
Unpert = Calc.TwoLevel(w0=w0,w=w,k_max=100,v=0,xi=1E-12);
Pert = Calc.TwoLevel(w0=w0,w=w,k_max=100,v=1E-4,xi=1E-2);
fine_tune(1) = struct(Vc=sqrt(w^2-del^2)+1E-15, delV=4E-4, nV=1001,...
    V_range=[]);
fine_tune(2) = struct(Vc=sqrt(4*w^2-del^2)+1E-15, delV=4E-4, nV=1001,...
    V_range=[]);
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
Details = Calc.Details(script='Calc_TwoLevel.m',...
    model='Two-level system',...
    calculation='Adiabatic Floquet eigenstate calculation',...
    objects=[Pert,Unpert],...
    variables=struct(V_range=V_range,fine_tune=fine_tune));
%% Calculate eigenstates
Res_Unpert_AE = Unpert.AEEigs(V_range);
Res_Unpert_QE = Unpert.QEEigs(V_range);
Res_Pert_AE = Pert.AEEigs(V_range);
Res_Pert_QE = Pert.QEEigs(V_range);
%% PostProcess
k_max2 = Unpert.k_max2;
for iV = 1:nV
    for iN = 1:2
        Res_Pert_AE(iV,iN).P = zeros(k_max2,1);
        Res_Pert_QE(iV,iN).P = zeros(k_max2,1);
        Res_Unpert_AE(iV,iN).P = zeros(k_max2,1);
        Res_Unpert_QE(iV,iN).P = zeros(k_max2,1);
        tPsi1 = reshape(Res_Pert_AE(iV,iN).Psi,2,[]);
        tPsi2 = reshape(Res_Pert_QE(iV,iN).Psi,2,[]);
        tPsi3 = reshape(Res_Unpert_AE(iV,iN).Psi,2,[]);
        tPsi4 = reshape(Res_Unpert_QE(iV,iN).Psi,2,[]);
        for ik = 1:k_max2
            Res_Pert_AE(iV,iN).P(ik) = tPsi1(:,ik)' * tPsi1(:,ik);
            Res_Pert_QE(iV,iN).P(ik) = tPsi2(:,ik)' * tPsi2(:,ik);
            Res_Unpert_AE(iV,iN).P(ik) = tPsi3(:,ik)' * tPsi3(:,ik);
            Res_Unpert_QE(iV,iN).P(ik) = tPsi4(:,ik)' * tPsi4(:,ik);
        end
    end
end
%% Save results
save('Calc_TwoLevel.mat',...
    'Res_Unpert_QE','Res_Unpert_AE',...
    'Res_Pert_QE','Res_Pert_AE',...
    'Details',...
    '-v7.3');