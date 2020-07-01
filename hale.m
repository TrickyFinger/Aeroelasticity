
L = mass * g; % lift force
EI = 2e4;

%% Mass parameters (per unit span)
m        = 0.75;
Icg_theta = 0.1;
Icg_beta  = .01;
W         = 30;

%% Geometric parameters
b       = 0.5;
a       = 0;
c       = 1.0;
x_theta = 0;
% x_beta  = 0;
S       = c*1;

%% Stiffness parameters
K_h     = 2818.8;
K_theta = 37.3;
K_beta  = 0;

%% Define output
parameters.m       = ma+mf;
% parameters.S_theta = ma*x_theta*b+mf*(c-a+x_beta)*b;
parameters.S_theta = 0.08587;
parameters.S_beta  = mf*x_beta*b;
% parameters.I_theta = ma*(x_theta*b)^2+mf*(c-a+x_beta)^2*b^2+Icg_theta+Icg_beta;
parameters.I_theta = 0.01347;
parameters.I_beta  = mf*(x_beta*b)^2+Icg_beta;
parameters.W       = W;
parameters.a       = a;
parameters.b       = b;
parameters.c       = c;
parameters.S       = S;
parameters.K_h     = K_h;
parameters.K_theta = K_theta; 
parameters.K_beta  = K_beta;
