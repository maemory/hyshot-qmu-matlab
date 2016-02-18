%% Import and Combine different data sets
clear all; clc;

nDim = 7;
% FFR 25, 35, and 40 only have 14 samples (rest of values are zeros)
% FFR 30 has 50 samples
% FFR 35 has 14 samples and the min & max active variable samples
nSamples = 50;

X_all    = zeros(nDim,nSamples);
Xhat_all = X_all;

% import intitial 14 samples
load first14_samples.mat X_l X_u X Xhat

X_all(:,1:14)    = X;
Xhat_all(:,1:14) = Xhat;

clear X Xhat

% import next 36 samples
load second36_samples.mat X Xhat

X_all(:,15:50)    = X;
Xhat_all(:,15:50) = Xhat;

clear X Xhat
X = X_all;
Xhat = Xhat_all;
clear X_all Xhat_all

load result_data.mat

%% Active Subspace Post Processing

% to match Paul's ordering of input variables
% our data:    Ptot Htot  AoA      turbI length_t t_ramp t_cowl
% Paul's data: AoA  turbI length_t Ptot  Htot     t_cowl t_ramp
indexed = [3, 4, 5, 1, 2, 7, 6];

FFR = 25;
metric = 'EP';

switch FFR
    case 25
        M = 14;
        
        switch metric
            case 'EP'
                QoI = EP_25;
            case 'SV'
                QoI = 100*ones(size(SV_25)) - SV_25;
            otherwise
                disp('unrecognized metric');
        end
        
        EquivRatio = ER_25;
        
    case 30
        M = 14;
        
        switch metric
            case 'EP'
                QoI = EP_30;
            case 'SV'
                QoI = 100*ones(size(SV_30)) - SV_30;
            otherwise
                disp('unrecognized metric');
        end
        
        EquivRatio = ER_30;

        
    case 35
        M = 14;
        
        switch metric
            case 'EP'
                QoI = EP_35;
            case 'SV'
                QoI = 100*ones(size(SV_35)) - SV_35;
            otherwise
                disp('unrecognized metric');
        end
        
        EquivRatio = ER_35;
        
    case 40
        M = 14;
        
        switch metric
            case 'EP'
                QoI = EP_40;
            case 'SV'
                QoI = 100*ones(size(SV_40)) - SV_40;
            otherwise
                disp('unrecognized metric');
        end
        
        EquivRatio = ER_40;
        
    otherwise
        disp('unrecognized FFR');
end

% compute weights (hard code that min/max generated from first 14 samples)
uhat = [ones(14,1) Xhat(:,1:14)'] \ QoI(1:14);
w = uhat(2:nDim+1) / norm(uhat(2:nDim+1));

disp(sprintf('FFR%d weights:',FFR));
disp(sprintf('AoA \t\t%f\nturbI \t\t%f\nlength_t \t%f\nPtot \t\t%f\nHtot \t\t%f\nt_cowl \t\t%f\nt_ramp \t\t%f',w(indexed)));

% compute min/max range values
A = [];
b = [];
Aeq = [];
beq = [];
lb = -1 * ones(nDim,1);
ub = 1* ones(nDim,1);
x0 = zeros(nDim,1);

x_min = fmincon(@(x) (w' * x),x0,A,b,Aeq,beq,lb,ub);
x_max = fmincon(@(x) -1*(w' * x),x0,A,b,Aeq,beq,lb,ub);
if (FFR == 30)
    Xhat_min = x_min;
    Xhat_max = x_max;
    X_min = 0.5 * ((X_u - X_l) .* x_min + (X_u + X_l));
    X_max = 0.5 * ((X_u - X_l) .* x_max + (X_u + X_l));
elseif (FFR == 35)
    Xhat_min = x_min;
    Xhat_max = x_max;
    X_min = 0.5 * ((X_u - X_l) .* x_min + (X_u + X_l));
    X_max = 0.5 * ((X_u - X_l) .* x_max + (X_u + X_l));
end
disp('Max active variable datum (Paul''s variable ordering)');
disp(sprintf('AoA \t\t%f\nturbI \t\t%f\nlength_t \t%f\nPtot \t\t%f\nHtot \t\t%f\nt_cowl \t\t%f\nt_ramp \t\t%f',x_max(indexed)));


% build SSPlot
figure(1)
plot(Xhat(:,1:M)' * w, QoI(1:M), 'o',...
    'MarkerSize', 10,...
    'LineWidth', 1,...
    'MarkerEdgeColor', [0.5, 0.5, 0.5],...
    'MarkerFaceColor', [0.5, 0.5, 0.5]);
if (FFR == 30)
   hold on
       plot([Xhat_min,Xhat_max]' * w, QoI(51:52), 'o',...
        'MarkerSize', 10,...
        'LineWidth', 1,...
        'MarkerEdgeColor', [0.9, 0.5, 0.5],...
        'MarkerFaceColor', [0.9, 0.5, 0.5]);
   hold off
elseif (FFR == 35)
   hold on
       plot([Xhat_min,Xhat_max]' * w, QoI(15:16), 'o',...
        'MarkerSize', 10,...
        'LineWidth', 1,...
        'MarkerEdgeColor', [0.9, 0.5, 0.5],...
        'MarkerFaceColor', [0.9, 0.5, 0.5]);
   hold off
end
ax = gca;
ax.FontSize = 18;
ax.FontName = 'Times New Roman';
ax.Box = 'off';
ax.TickDir = 'out';

figure(2)
plot(Xhat(:,1:M)' * w, EquivRatio(1:M), 'o',...
    'MarkerSize', 10,...
    'LineWidth', 1,...
    'MarkerEdgeColor', [0.5, 0.5, 0.5],...
    'MarkerFaceColor', [0.5, 0.5, 0.5]);
ax = gca;
ax.FontSize = 18;
ax.FontName = 'Times New Roman';
ax.Box = 'off';
ax.TickDir = 'out';

%% Min/Max AS samples

if (0)
    % compute min/max range values
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = -1 * ones(nDim,1);
    ub = 1* ones(nDim,1);
    x0 = zeros(nDim,1);

    x_min = fmincon(@(x) (w' * x),x0,A,b,Aeq,beq,lb,ub);
    av_min = x_min' * w;
    disp(sprintf('Min active variable datum: %f',av_min));
    disp(sprintf('Ptot \t\t%d\nHtot \t\t%d\nAoA \t\t%d\nturbI \t\t%d\nlength_t \t%d\nt_ramp \t\t%d\nt_cowl \t\t%d',round(x_min)));


    x_max = fmincon(@(x) -1*(w' * x),x0,A,b,Aeq,beq,lb,ub);
    av_max = x_max' * w;
    disp(sprintf('Max active variable datum: %f',av_max));
    disp(sprintf('Ptot \t\t%d\nHtot \t\t%d\nAoA \t\t%d\nturbI \t\t%d\nlength_t \t%d\nt_ramp \t\t%d\nt_cowl \t\t%d',round(x_max)));

    % 25 max 1.975614 [1 -1 1 1 1 1 -1]
    % 30 max 1.995537 [1 -1 1 1 -1 1 1]
    % 35 max 1.993194 [1 -1 1 1 1 1 -1]
    % 40 max 1.773432 [-1 -1 1 1 -1 1 1]

    %% Compute BCs for 2D min/max sims
    X = 0.5 * ((X_u - X_l) .* x_min + (X_u + X_l));
    X_min = X;
    Xhat_min = x_min;
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

    t_ramp = X(6,:);
    t_cowl = X(7,:) + 0.291;

    disp(' ');
    disp('Min parameters');
    disp('--------------');
    disp(sprintf('Mach : %f',Mach));
    disp(sprintf('AoA  : %f',AoA));
    disp(sprintf('Pstat: %f',Pstat));
    disp(sprintf('Tstat: %f',Tstat));
    disp(sprintf('rho  : %f',rho_stat));
    disp(sprintf('Ux   : %f',Ux));
    disp(sprintf('Uy   : %f',Uy));
    disp(sprintf('k    : %f',k));
    disp(sprintf('omega: %f',omega));
    disp(sprintf('t_cwl: %f',t_cowl));
    disp(sprintf('t_rmp: %f',t_ramp));

    X = 0.5 * ((X_u - X_l) .* x_max + (X_u + X_l));
    X_max = X;
    Xhat_max = x_max;
    
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

    t_ramp = X(6,:);
    t_cowl = X(7,:) + 0.291;
    disp(' ');
    disp('Max parameters');
    disp('--------------');
    disp(sprintf('Mach : %f',Mach));
    disp(sprintf('AoA  : %f',AoA));
    disp(sprintf('Pstat: %f',Pstat));
    disp(sprintf('Tstat: %f',Tstat));
    disp(sprintf('rho  : %f',rho_stat));
    disp(sprintf('Ux   : %f',Ux));
    disp(sprintf('Uy   : %f',Uy));
    disp(sprintf('k    : %f',k));
    disp(sprintf('omega: %f',omega));
    disp(sprintf('t_cwl: %f',t_cowl));
    disp(sprintf('t_rmp: %f',t_ramp));
end