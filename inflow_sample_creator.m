%% Active Subspace routine for HyShot II QMU studies
clear all; clc;

% dimensions
m = 7;
% parameter ranges
% m=7: Ptot Htot AoA turbI length_t t_ramp t_cowl
X_l = [16.448 3.0551 2.6 0.001 0.1325 0.087 0.030]';
X_n = [17.730 3.2415 3.6 0.01 0.245 0.145 0.050]';
X_u = [19.012 3.428 4.6 0.019 0.3575 0.203 0.070]';

% number of samples
M = 36;

% generate sample values for input parameters
Xhat = 2*rand(m,M) - 1;
X = 0.5*(repmat(X_u - X_l,1,M).*Xhat + repmat(X_u + X_l,1,M));


%% Format for use in Joe.in

%X = X_n;

% fluid properties:
gamma = 1.4;
R = 287;

Pstat = X(1,:) * 1e6 * 1.16e-4;

Ttot = X(2,:) * 1e6 * 6.8718e-4 + 508.136;
Tstat = Ttot * 0.0978;

rho_stat = Pstat./(R * Tstat);

AoA = -X(3,:);

Umag = 1.332 * sqrt(X(2,:) * 1e6);
Ux = cosd(AoA) .* Umag;
Uy = sind(AoA) .* Umag;

sos = sqrt(gamma * R * Tstat);
Mach = Umag ./ sos;

k = 1.5 * (Umag .* X(4,:)).^2;

omega = sqrt(k) ./ ((0.09^0.25) * X(5,:));

t_ramp = X(6,:)
t_cowl = X(7,:) + 0.291



%% Equivalence ratio
% equation is linear regression from DLR data
ER = [0.25 0.3 0.35 0.4]
P_injector = 1606916.78 * ER



