clc; clear all;
w=8.3;
xi=1E-2; tol=1E-12;
obj = Calc.BoxParticle(20,w=w,k_max=200,xi=1E-12);
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
Details = Calc.Details(script='Calc_BoxParticle.m',...
    model='Particle in a box',...
    calculation='Adiabatic Floquet eigenstate calculation',...
    objects=obj,...
    variables=struct(V_range=V_range,fine_tune=fine_tune));
%% Calculate eigenstates
Res_AE = obj.AEEigs(V_range);
Res_QE = obj.QEEigs(V_range);
%% PostProcess
k_max2 = obj.k_max2;
for iV = 1:nV
    for iN = 1:2
        Res_AE(iV,iN).P = zeros(k_max2,1);
        Res_QE(iV,iN).P = zeros(k_max2,1);
        tPsi1 = reshape(Res_AE(iV,iN).Psi,obj.N,[]);
        tPsi2 = reshape(Res_QE(iV,iN).Psi,obj.N,[]);
        for ik = 1:k_max2
            Res_AE(iV,iN).P(ik) = tPsi1(:,ik)' * tPsi1(:,ik);
            Res_QE(iV,iN).P(ik) = tPsi2(:,ik)' * tPsi2(:,ik);
        end
    end
end
%% Save results
save('Calc_BoxParticle.mat',...
    'Res_QE','Res_AE',...
    'Details',...
    '-v7.3');