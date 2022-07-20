clc; clear all;
close all;
%% Load data
Calc1=load('Calc_BoxParticle.mat');
Calc2=load('Calc_BoxParticle2.mat');
% w = Calc1.Details.objects(1).w;
% k_max2 = Calc1.Details.objects(1).k_max2;
% V_range = Calc1.Details.variables.V_range;
% nV = length(V_range);
%% Plot Calc2
plotObj = Plot.BoxParticle(Calc2.Details.objects(1));
[fg(1),ax(1)] = plotObj.Variational_Overlap(Calc2.Res_AE,1,...
    File='PlotB29.png',normalized=true);
[fg(2),ax(2)] = plotObj.Variational_varEps(Calc2.Res_AE,...
    File='PlotB29-1.png',exact=true);