clc; clear all;
close all;
%% Load data
load('CalcRes.mat');
%% Recreate missing object
obj = Calc.HubardDimerOrbital(t=1,v0=2,v1=5,xi=1E-4);
obj_ex = Calc.HubardDimerExact(t=1,v0=2,v1=5,xi=1E-4);
%%
plotObj = Plot.HubardDimer(obj);
[fg(1),ax(1)]=plotObj.Theta(Exact=res_Exact,FHF=res_FHF,NSCF=res_NSCF,...
    File='PlotB30.png');
[fg(2),ax(2)]=plotObj.AverageEnergy(Exact=res_Exact,FHF=res_FHF,...
    File='PlotB30-1.png');
[fg(3:5),ax(3:5)]=plotObj.OverlapEigen(NSCF=res_NSCF,FHF=res_FHF,merged=false,...
    Files={'PlotB31.png','PlotB31-1.png','PlotB31-2.png'});
%% Manual plots
[fg(6),ax(6)]=plotObj.prePlot;
plot_data=abs(res_tdHF.Res(16,1).S_th2);
plot_data=min(plot_data,1);
plot(ax(6),res_tdHF.t_range,plot_data,...
    LineWidth=1.0);
ylim(ax(6),[0 1+1E-3]);
xlabel(ax(6),'Time',...
    Interpreter='latex');
ylabel(ax(6),'Overlap $\vert\langle\Psi_{FHF,0}(t)\vert\Psi_{td-HF}(t)\rangle\vert$',...
    Interpreter='latex')
ax(6).XScale="log";
plotObj.postPlot(fg(6),ax(6),...
    File='PlotB32.png');
%%
[fg(7),ax(7)]=plotObj.prePlot;
plot_data=abs(res_tdHF.Res(1,1).S_th2);
% plot_data=min(plot_data,1);
plot(ax(7),res_tdHF.t_range,plot_data,...
    LineWidth=1.0);
% yl=
% ylim(ax(7),[ 1+1E-3]);
xlabel(ax(7),'Time',...
    Interpreter='latex');
ylabel(ax(7),'Overlap $\vert\langle\Psi_{FHF,0}(t)\vert\Psi_{td-HF}(t)\rangle\vert$',...
    Interpreter='latex')
ax(7).XScale="log";
plotObj.postPlot(fg(7),ax(7),...
    File='PlotB32-1.png');