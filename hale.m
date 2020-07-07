clear all;

target_freq_bending = 2.247;
target_freq_torsion = 31.146;

%% Tune the location of the typical section
weight_wb = 0.8;
weight_wt = 1 - weight_wb;
min_TOL = 1;
for ts_frac = 0.5:0.001:0.9
    TS = getTypicalSectionParam(ts_frac);
    [Ms, Ks] = getTypicalSectionMatrices(TS);
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
% reinitialize the TS parameters
TS = getTypicalSectionParam(loc_ts);
[Ms, Ks] = getTypicalSectionMatrices(TS);
[V, D] = eig(Ks, Ms);
eig_freq = sqrt(D);
freq_bending = eig_freq(1);
freq_torsion = eig_freq(4);

%% Static model divergence speed
% flow = getFlowParameters;
% q_div = TS.K_theta / (TS.e * TS.c * flow.C_L_alpha * TS.S);
% speed_div = sqrt(2*q_div/flow.rho);
% 
% i = 1;
% for q = 0:0.5:50
%     a0 = TS.K_h * (TS.K_theta - TS.e*TS.c*q*TS.S*flow.C_L_alpha);
%     a2 = TS.m * TS.K_theta + TS.I_theta * TS.K_h ... 
%          - (2*TS.m*TS.e*TS.c+TS.S_theta) * q * TS.S * flow.C_L_alpha;
%     a4 = TS.m * TS.I_theta - TS.S_theta^2;
%     p(i, :) = roots([a4 0 a2 0 a0]);
%     for j = 1:4
%         if imag(p(i, j)) >= 0
%             p_real(i, j) = real(p(i, j));
%             p_imag(i, j) = imag(p(i, j));
%         else
%             p_real(i, j) = 0;
%             p_imag(i, j) = 0;
%         end
%     end
%     i = i + 1;
% end
% figure(1)
% hold on
% plot((0:0.5:50), p_real(:, 1));
% plot((0:0.5:50), p_real(:, 2));
% plot((0:0.5:50), p_real(:, 3));
% plot((0:0.5:50), p_real(:, 4));
% hold off
% figure(2)
% hold on
% plot((0:0.5:50), p_imag(:, 1));
% plot((0:0.5:50), p_imag(:, 2));
% plot((0:0.5:50), p_imag(:, 3));
% plot((0:0.5:50), p_imag(:, 4));
% hold off

%% Quasi-steady model flutter
% a0 = TS.K_h * (TS.K_theta - TS.e*TS.c*q*TS.S*flow.C_L_alpha);
% a1 = 
% a2 = TS.m * TS.K_theta + TS.I_theta * TS.K_h ...
%     - (2*TS.m*TS.e*TS.c+TS.S_theta) * q * TS.S * flow.C_L_alpha;
% a4 = TS.m * TS.I_theta - TS.S_theta^2;
% p(i, :) = roots([a4 0 a2 0 a0]);
% Aerodynamic stiffness matrix
Ka = [0                 -flow.C_L_alpha;
      0    (0.5 + TS.a) * TS.b * TS.S * flow.C_L_alpha];
i = 1;
for q = 0:0.5:50
    v = sqrt(2*q/flow.rho);
    
    % Aerodynamic damping matrix
    Ca = (q/v).*[-TS.S * flow.C_L_alpha             -TS.S * TS.b * (pi+2*pi*(0.5-TS.a));
              0.5 * TS.b * TS.S * flow.C_L_alpha                   0                ];
    p = sym('p');
    poly = det(Ms.*(p*p) - Ca.*p + (Ks-Ka));
    coeffs = double(coeffs(poly));
    a4 = coeffs(1);
    a3 = coeffs(2);
    a2 = coeffs(3);
    a1 = coeffs(4);
    a0 = coeffs(5);
    p(i, :) = roots([a4 a3 a2 a1 a0]);
    for j = 1:4
        if imag(p(i, j)) >= 0
            p_real(i, j) = real(p(i, j));
            p_imag(i, j) = imag(p(i, j));
        else
            p_real(i, j) = 0;
            p_imag(i, j) = 0;
        end
    end
    i = i + 1;
end
figure(3)
hold on
plot((0:0.5:50), p_real(:, 1));
plot((0:0.5:50), p_real(:, 2));
plot((0:0.5:50), p_real(:, 3));
plot((0:0.5:50), p_real(:, 4));
hold off
figure(4)
hold on
plot((0:0.5:50), p_imag(:, 1));
plot((0:0.5:50), p_imag(:, 2));
plot((0:0.5:50), p_imag(:, 3));
plot((0:0.5:50), p_imag(:, 4));
hold off

function TS = getTypicalSectionParam(ts_frac)
%% Mass parameters (per unit span)
m         = 0.75;
Icg_theta = 0.1;
W         = 30;

%% Geometric parameters
b       = 0.5;
a       = 0;
c       = 2 * b;
e       = 0.25;   % distance between AC and EA
x_theta = 0;
S       = c;
semi_span = 16;

%% Location of the typical section
loc_ts = ts_frac * semi_span;

%% Stiffness parameters
EI      = 2e4;
GJ      = 1e4;
K_h     = 3 * EI / loc_ts^3;
K_theta = GJ / loc_ts;

%% Define output
TS.m       = m;
TS.S_theta = m * x_theta * b;
TS.I_theta = m * (x_theta*b)^2 + Icg_theta;
TS.W       = W;
TS.a       = a;
TS.b       = b;
TS.c       = c;
TS.e       = e;
TS.S       = S;
TS.K_h     = K_h / loc_ts;
TS.K_theta = K_theta / loc_ts; 
end

function [Ms, Ks] = getTypicalSectionMatrices(param)
%% Typical section structural matrix
Ms = [param.m          param.S_theta;                                                     
      param.S_theta    param.I_theta];                                       
      

Ks = [param.K_h         0        ;
        0           param.K_theta];
end

function flow = getFlowParameters
rho       = 0.0889;
V         = 23;
% alpha0    = 5*pi/180;
C_M_AC    = 0;
C_L_alpha = 2*pi;

q = .5 * rho * V^2;

flow.rho       = rho;
flow.V         = V;
% flow.alpha0    = alpha0;
flow.C_M_AC    = C_M_AC;
flow.C_L_alpha = C_L_alpha;
flow.q         = q;
end