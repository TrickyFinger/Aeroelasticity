clear all;

target_freq_bending = 2.247;
target_freq_torsion = 31.146;

%% Tune the location of the typical section
weight_wb = 0.8;
weight_wt = 1 - weight_wb;
min_TOL = 1;
for ts_frac = 0.5:0.001:0.9
    param = getTypicalSectionParam(ts_frac);
    [Ms, Ks] = getTypicalSectionMatrices(param);
    [V, D] = eig(Ks, Ms);
    eig_freq = sqrt(D);
    freq_bending = eig_freq(1);
    freq_torsion = eig_freq(4);
    err_wb = abs((freq_bending - target_freq_bending) / target_freq_bending);
    err_wt = abs((freq_torsion - target_freq_torsion) / target_freq_torsion);
    TOL = weight_wb * err_wb + weight_wt * err_wt;
    if TOL < min_TOL
        min_TOL = TOL;
        loc_ts = ts_frac;
    end
end
param = getTypicalSectionParam(loc_ts);
[Ms, Ks] = getTypicalSectionMatrices(param);
[V, D] = eig(Ks, Ms);
eig_freq = sqrt(D);
freq_bending = eig_freq(1);
freq_torsion = eig_freq(4);



function param = getTypicalSectionParam(ts_frac)
%% Mass parameters (per unit span)
m        = 0.75;
Icg_theta = 0.1;
W         = 30;

%% Geometric parameters
b       = 0.5;
a       = 0;
c       = 2 * b;
x_theta = 0;
S       = c;
semi_span = 16;

%% Location of the typical section
loc_ts = ts_frac * semi_span;

%% Stiffness parameters
EI = 2e4;
GJ = 1e4;
K_h     = 3 * EI / loc_ts^3;
K_theta = GJ / loc_ts;

%% Define output
param.m       = m;
param.S_theta = m * x_theta * b;
param.I_theta = m * (x_theta*b)^2 + Icg_theta;
param.W       = W;
param.a       = a;
param.b       = b;
param.c       = c;
param.S       = S;
param.K_h     = K_h / loc_ts;
param.K_theta = K_theta / loc_ts; 

end

function [Ms, Ks] = getTypicalSectionMatrices(param)
%% Typical section structural matrix
Ms = [param.m          param.S_theta;                                                     
      param.S_theta    param.I_theta];                                       
      

Ks = [param.K_h         0        ;
        0           param.K_theta];
end