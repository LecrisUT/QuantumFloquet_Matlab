clc; clear all;
close all;
load('CalcRes.mat');
plotObj = Plot.HubardDimer;
plotObj.Theta(Exact=res_Exact,FHF=res_FHF,NSCF=res_NSCF);
plotObj.AverageEnergy(Exact=res_Exact,FHF=res_FHF);
plotObj.OverlapEigen('NSCF',res_NSCF,'FHF',res_FHF,'merged',false);
plotObj.OverlapEigen('NSCF',res_NSCF,'FHF',res_FHF,'merged',true);
plotObj.OverlapProp(1:4,NSCF=res_NSCF,FHF=res_FHF,...
    plot_S_th=true,plot_S_min=false);
plotObj.OverlapProp(1:4,NSCF=res_NSCF,FHF=res_FHF,...
    plot_S_th=false,plot_S_min=true);
% plotObj.OverlapProp(4,NSCF=res_NSCF,FHF=res_FHF,...
%     plot_S_th=true,plot_S_min=false);
% plotObj.OverlapProp(4,NSCF=res_NSCF,FHF=res_FHF,...
%     plot_S_th=false,plot_S_min=true);
plotObj.t_max(4,'NSCF',res_NSCF,'FHF',res_FHF,'merged',false);
%%
import mlreportgen.report.*
import mlreportgen.dom.*
rpt = Report('Hubbard_Dimer','pdf');
for curr_fig = plotObj.figs
    fig = Figure(curr_fig);
    rpt.append(fig);
end
close(rpt);
rptview(rpt);