% active subspace post-processing
clear all; clc;

load /Users/memory/Dropbox/Matlab/Hyshot_QMU/first14_samples.mat

ER = zeros(14,4);
EP = zeros(14,4);
SV = zeros(14,4);

% ER = 0.25
EP(:,1) = [
2.20546E+00
2.03850E+00
2.06730E+00
2.10156E+00
2.10488E+00
2.06879E+00
2.23954E+00
2.04173E+00
2.15578E+00
2.12349E+00
2.10296E+00
2.25631E+00
2.17452E+00
2.16016E+00
];

SV(:,1) = [
9.90434E+01
9.89469E+01
9.90757E+01
9.90607E+01
9.89707E+01
9.89782E+01
9.90637E+01
9.90214E+01
9.89135E+01
9.91367E+01
9.90012E+01
9.89550E+01
9.88464E+01
9.90519E+01
];

% ER = 0.30
ER(:,2) = [
2.8050E-01
3.0966E-01
2.9883E-01
3.0320E-01
3.0468E-01
2.9506E-01
2.7979E-01
3.1435E-01
3.1678E-01
2.9634E-01
3.0462E-01
2.8189E-01
2.9934E-01
2.8789E-01
];

EP(:,2) = [
2.3922E+00
2.2547E+00
2.2649E+00
2.3032E+00
2.3176E+00
2.2747E+00
2.4381E+00
2.2485E+00
2.4245E+00
2.3206E+00
2.3099E+00
2.4574E+00
2.4247E+00
2.3524E+00
];

SV(:,2) = [
9.9003E+01
9.8843E+01
9.9026E+01
9.8978E+01
9.8868E+01
9.8908E+01
9.8951E+01
9.8935E+01
9.8722E+01
9.9044E+01
9.8907E+01
9.8908E+01
9.8728E+01
9.8997E+01
];

% ER = 0.35

EP(:,3) = [
2.6267E+00
2.5484E+00
2.4837E+00
2.5603E+00
2.6084E+00
2.5283E+00
2.6950E+00
2.5320E+00
2.7107E+00
2.5463E+00
2.5956E+00
2.7359E+00
2.7028E+00
2.5888E+00
];

SV(:,3) = [
9.8866E+01
9.8522E+01
9.8841E+01
9.8798E+01
9.8588E+01
9.8642E+01
9.8736E+01
9.8704E+01
9.8289E+01
9.8847E+01
9.8648E+01
9.8725E+01
9.8333E+01
9.8821E+01
];

% ER = 0.40

EP(:,4) = [
2.9092E+00
2.9167E+00
2.7792E+00
2.8325E+00
2.9505E+00
2.8252E+00
2.9693E+00
2.8589E+00
3.1294E+00
2.8231E+00
2.9120E+00
3.0410E+00
3.1373E+00
2.8712E+00
];

SV(:,4) = [
9.8581E+01
9.7451E+01
9.8480E+01
9.8295E+01
9.7809E+01
9.7998E+01
9.8354E+01
9.8009E+01
9.7040E+01
9.8514E+01
9.8010E+01
9.8272E+01
9.7079E+01
9.8544E+01
];

%% Load Paul's data

data_in = importdata('./HyShotII.txt');

% place Pauls ordering in ours (for min max)
% Ptot Htot AoA turbI length_t t_ramp t_cowl
rev_indexed = [4, 5, 1, 2, 3, 7, 6];

% to match Paul's ordering of input variables
% AoA turbI length_t Ptot Htot t_cowl t_ramp
indexed = [3, 4, 5, 1, 2, 7, 6];

XP = data_in.data(:,2:8);
XP(:,4:5) = XP(:,4:5) .* 1e-6;
XP(:,2) = XP(:,2) .* 1e-2;
XP(:,3) = XP(:,3) .* 1e-3;
XP = XP(:,rev_indexed)';

XPhat = (2 * XP - repmat(X_u + X_l,1,M)) ./ repmat(X_u - X_l,1,M);

EP_P = zeros(M,4);
w_ep_P = zeros(m,4);

EP_P(:,2) = data_in.data(:,9);
EP_P(:,3) = data_in.data(:,10);

% compute weights for Paul's data
uhat = [ones(M,1) XPhat'] \ EP_P(:,2);
w_ep_P(:,2) = uhat(2:m+1) / norm(uhat(2:m+1));

uhat = [ones(M,1) XPhat'] \ EP_P(:,3);
w_ep_P(:,3) = uhat(2:m+1) / norm(uhat(2:m+1));

%% compute linear regression weights
close all;

for FFR = 1:4
    uhat = [ones(M,1) Xhat'] \ EP(:,FFR);
    w_ep(:,FFR) = uhat(2:m+1) / norm(uhat(2:m+1));
    
    subV(:,FFR) = 100*ones(size(SV(:,FFR))) - SV(:,FFR);
    uhat = [ones(M,1) Xhat'] \ subV(:,FFR);
    w_sv(:,FFR) = uhat(2:m+1) / norm(uhat(2:m+1));

    % same weights ordering as Paul's paper
    disp(sprintf('FFR%d w exit press:',FFR*5+20))
    disp(w_ep(indexed,FFR))

    % same weights ordering as Paul's paper
    disp(sprintf('FFR%d w subsonic v:',FFR*5+20))
    disp(w_sv(indexed,FFR))
    
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = -1 * ones(m,1);
    ub = 1* ones(m,1);
    x0 = zeros(m,1);

%     x_min = fmincon(@(x) w_ep(:,FFR)' * x,x0,A,b,Aeq,beq,lb,ub);
%     av_min = x_min' * w_ep(:,FFR);
%     disp('EP min a.v.:')
%     disp(av_min)
%     disp('EP min weights:')
%     disp(x_min')

    x_max = fmincon(@(x) -1*(w_ep(:,FFR)' * x),x0,A,b,Aeq,beq,lb,ub);
    av_max = x_max' * w_ep(:,FFR);

    disp(sprintf('FFR%d max(AV):',FFR*5+20))
    disp(av_max)
    disp(sprintf('FFR%d max sample:',FFR*5+20))
    disp(x_max')
end

% for plotting
FFR = 1;

% sufficient summary plot (exit pressure)
figure(1)
plot(Xhat' * w_ep(:,FFR), EP(:,FFR), 'o',...
    'MarkerSize', 10,...
    'LineWidth', 1,...
    'MarkerEdgeColor', [0.5, 0.5, 0.5],...
    'MarkerFaceColor', [0.5, 0.5, 0.5]);
if (FFR == 2 || FFR == 3)
    hold on
    plot(XPhat' * w_ep_P(:,FFR), EP_P(:,FFR), 'o',...
        'MarkerSize', 10,...
        'LineWidth', 1,...
        'MarkerEdgeColor', [24,87,155]./255,...
        'MarkerFaceColor', [24,87,155]./255);
    hold off
end
% axis([-0.9 0.5 2.2 2.5]);
ax = gca;
ax.FontSize = 18;
ax.FontName = 'Times New Roman';
ax.Box = 'off';
ax.TickDir = 'out';

% sufficient summary plot (supersonic volume)
figure(2)
plot(Xhat' * w_sv(:,FFR), subV(:,FFR), 'o',...
    'MarkerSize', 10,...
    'LineWidth', 1,...
    'MarkerEdgeColor', [0.5, 0.5, 0.5],...
    'MarkerFaceColor', [0.5, 0.5, 0.5]);
% axis([-0.7 0.6 0.9 1.3]);
ax = gca;
ax.FontSize = 18;
ax.FontName = 'Times New Roman';
ax.Box = 'off';
ax.TickDir = 'out';


% bar chart evaluating weights
figure(3)
b = bar(cat(2,w_ep(indexed,FFR),w_ep_P(indexed,FFR)),1,...
            'EdgeColor', 'w');
        
b(1).FaceColor = [0.5, 0.5, 0.5];
b(2).FaceColor = [24,87,155]./255;
ax = gca;
ax.FontSize = 20;
ax.FontName = 'Times New Roman';
ax.Box = 'off';
ax.TickDir = 'out';
ax.XTickLabel = {'\alpha', 'I', 'L_t', 'P_0', 'H_0', 'x_{t,C}', 'x_{t,R}'};
ax.XTickLabelRotation = 0;
fig3 = gcf;
fig3.OuterPosition(3) = fig3.OuterPosition(3) + 50;

% compare between FFR
figure(4)
b = bar(cat(2,w_ep(indexed,1),w_ep(indexed,2),w_ep(indexed,3),w_ep(indexed,4)),1,...
            'EdgeColor', 'w');
b(1).FaceColor = 'k';
b(2).FaceColor = [149,149,149]./255;
b(3).FaceColor = [24,87,155]./255;
b(4).FaceColor = [134,188,247]./255;

ax = gca;
ax.FontSize = 20;
ax.FontName = 'Times New Roman';
ax.Box = 'off';
ax.TickDir = 'out';
ax.XTickLabel = {'\alpha', 'I', 'L_t', 'P_0', 'H_0', 'x_{t,C}', 'x_{t,R}'};
ax.XTickLabelRotation = 0;
fig3 = gcf;
fig3.OuterPosition(3) = fig3.OuterPosition(3) + 50;

% figure(4)
% plot(Xhat' * w_ep,ER30,'o',...
%     'MarkerSize', 10,...
%     'LineWidth', 1,...
%     'MarkerEdgeColor', 'k',...
%     'MarkerFaceColor', [0.5, 0.5, 0.5]);

%%

X = 0.5* ( (X_u - X_l).*[-1,-1,1,1,-1,1,1]' + (X_u + X_l));

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

Mach
AoA
Pstat
Tstat
rho_stat
Ux
Uy
k
omega
t_cowl
t_ramp


