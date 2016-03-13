% active subspace post-processing
clear all; clc;

M = 24;
m = 12;

MAX_HR = zeros(M+4,4);
INT_HR = zeros(M+4,4);

MAX_HR(:,1) = [
2.09543E+11
2.18537E+11
2.25695E+11
2.47187E+11
2.03189E+11
2.39773E+11
2.09532E+11
1.95636E+11
2.47692E+11
2.36010E+11
2.01018E+11
2.46948E+11
2.21769E+11
2.29203E+11
2.36624E+11
2.35324E+11
2.34909E+11
2.37980E+11
2.14693E+11
2.23315E+11
2.04919E+11
2.45784E+11
2.51214E+11
2.39258E+11
1.74881E+11
2.78825E+11
1.76175E+11
2.76149E+11
];

MAX_HR(:,2) = [
2.12910E+11
2.29582E+11
2.26642E+11
2.51225E+11
2.13684E+11
2.45333E+11
2.18642E+11
2.06436E+11
2.54975E+11
2.40264E+11
2.07980E+11
2.51943E+11
2.23920E+11
2.32731E+11
2.46812E+11
2.44077E+11
2.44031E+11
2.45696E+11
2.18786E+11
2.27954E+11
2.13160E+11
2.51083E+11
2.60306E+11
2.43964E+11
1.84369E+11
2.90364E+11
1.86003E+11
2.87873E+11
];

MAX_HR(:,3) = [
2.15737E+11
2.35519E+11
2.28821E+11
2.52654E+11
2.20146E+11
2.45737E+11
2.25573E+11
2.12363E+11
2.63673E+11
2.44073E+11
2.13464E+11
2.53492E+11
2.28177E+11
2.36594E+11
2.51590E+11
2.49405E+11
2.52038E+11
2.54208E+11
2.22196E+11
2.32320E+11
2.20232E+11
2.52165E+11
2.67270E+11
2.44563E+11
1.83378E+11
2.96513E+11
1.84898E+11
2.94027E+11
];

MAX_HR(:,4) = [
2.12836E+11
2.28429E+11
2.26636E+11
2.52700E+11
2.13104E+11
2.45872E+11
2.19754E+11
2.05002E+11
2.58548E+11
2.40482E+11
2.07509E+11
2.54287E+11
2.25475E+11
2.32940E+11
2.47386E+11
2.45138E+11
2.45242E+11
2.48241E+11
2.20043E+11
2.28447E+11
2.14554E+11
2.52658E+11
2.63216E+11
2.44534E+11
1.82081E+11
2.91837E+11
1.83801E+11
2.89853E+11
];

INT_HR(:,1) = [
1.54148E+05
1.66468E+05
1.68528E+05
1.63965E+05
1.57435E+05
1.51547E+05
1.63405E+05
1.58034E+05
1.72006E+05
1.57523E+05
1.47609E+05
1.66012E+05
1.59645E+05
1.52564E+05
1.78583E+05
1.74085E+05
1.70637E+05
1.64753E+05
1.65835E+05
1.68006E+05
1.57026E+05
1.56624E+05
1.81974E+05
1.66409E+05
1.12562E+05
1.72764E+05
1.14286E+05
1.72285E+05
];

INT_HR(:,2) = [
1.82265E+05
1.96857E+05
1.99086E+05
1.93967E+05
1.85925E+05
1.79938E+05
1.92477E+05
1.86667E+05
2.03428E+05
1.86550E+05
1.74947E+05
1.96398E+05
1.88741E+05
1.80922E+05
2.10981E+05
2.05824E+05
2.02027E+05
1.95067E+05
1.95776E+05
1.98863E+05
1.84832E+05
1.85750E+05
2.14886E+05
1.97174E+05
1.33148E+05
2.04301E+05
1.35118E+05
2.03707E+05
];

INT_HR(:,3) = [
1.98618E+05
2.14574E+05
2.16575E+05
2.11025E+05
2.02915E+05
1.96457E+05
2.09589E+05
2.03864E+05
2.21153E+05
2.03215E+05
1.91133E+05
2.13670E+05
2.05549E+05
1.97338E+05
2.29205E+05
2.23770E+05
2.20012E+05
2.12404E+05
2.13136E+05
2.16544E+05
2.01387E+05
2.02580E+05
2.33287E+05
2.14793E+05
1.42678E+05
2.18012E+05
1.44862E+05
2.17373E+05
];

INT_HR(:,4) = [
2.19898E+05
2.36458E+05
2.38662E+05
2.33334E+05
2.23758E+05
2.18435E+05
2.30325E+05
2.24344E+05
2.44065E+05
2.25332E+05
2.12231E+05
2.36156E+05
2.27464E+05
2.19300E+05
2.51948E+05
2.46441E+05
2.42721E+05
2.35002E+05
2.34985E+05
2.38802E+05
2.20851E+05
2.24625E+05
2.56217E+05
2.37321E+05
1.55203E+05
2.37245E+05
1.58178E+05
2.36533E+05
];

%% Normalized Data
% Xhat: samples normalized by min/max values [-1,1]

% X: raw sample values
% X_u: vector of max sample values (upper)
% X_l: vector of min sample values (lower)
load /Users/memory/Dropbox/Matlab/Hyshot_QMU/PAUL-chem24_samples.mat

Xhat = (2 * X - repmat(X_u + X_l,1,M)) ./ repmat(X_u - X_l,1,M);

%% compute linear regression weights
close all;

% FFR25
uhat = [ones(M,1) Xhat'] \ MAX_HR(1:M,1);
w_maxHR_25 = uhat(2:m+1) / norm(uhat(2:m+1));

uhat = [ones(M,1) Xhat'] \ INT_HR(1:M,1);
w_intHR_25 = uhat(2:m+1) / norm(uhat(2:m+1));

% FFR30
uhat = [ones(M,1) Xhat'] \ MAX_HR(1:M,2);
w_maxHR_30 = uhat(2:m+1) / norm(uhat(2:m+1));

uhat = [ones(M,1) Xhat'] \ INT_HR(1:M,2);
w_intHR_30 = uhat(2:m+1) / norm(uhat(2:m+1));

% FFR35
uhat = [ones(M,1) Xhat'] \ MAX_HR(1:M,3);
w_maxHR_35 = uhat(2:m+1) / norm(uhat(2:m+1));

uhat = [ones(M,1) Xhat'] \ INT_HR(1:M,3);
w_intHR_35 = uhat(2:m+1) / norm(uhat(2:m+1));

% FFR40
uhat = [ones(M,1) Xhat'] \ MAX_HR(1:M,4);
w_maxHR_40 = uhat(2:m+1) / norm(uhat(2:m+1));

uhat = [ones(M,1) Xhat'] \ INT_HR(1:M,4);
w_intHR_40 = uhat(2:m+1) / norm(uhat(2:m+1));

%% minimization
    
A = [];
b = [];
Aeq = [];
beq = [];
lb = -1 * ones(m,1);
ub = 1* ones(m,1);
x0 = zeros(m,1);

% min MHR
xhat_min_MHR_25 = fmincon(@(x) 1*(w_maxHR_25' * x),x0,A,b,Aeq,beq,lb,ub);
av_min = xhat_min_MHR_25' * w_maxHR_25;
Xmin_MHR_25 = 0.5 * (xhat_min_MHR_25 .* (X_u - X_l) + (X_u + X_l));

xhat_min_MHR_30 = fmincon(@(x) 1*(w_maxHR_30' * x),x0,A,b,Aeq,beq,lb,ub);
av_min = xhat_min_MHR_30' * w_maxHR_30;
Xmin_MHR_30 = 0.5 * (xhat_min_MHR_30 .* (X_u - X_l) + (X_u + X_l));

xhat_min_MHR_35 = fmincon(@(x) 1*(w_maxHR_35' * x),x0,A,b,Aeq,beq,lb,ub);
av_min = xhat_min_MHR_35' * w_maxHR_35;
Xmin_MHR_35 = 0.5 * (xhat_min_MHR_35 .* (X_u - X_l) + (X_u + X_l));

xhat_min_MHR_40 = fmincon(@(x) 1*(w_maxHR_40' * x),x0,A,b,Aeq,beq,lb,ub);
av_min = xhat_min_MHR_40' * w_maxHR_40;
Xmin_MHR_40 = 0.5 * (xhat_min_MHR_40 .* (X_u - X_l) + (X_u + X_l));

% max MHR
xhat_max_MHR_25 = fmincon(@(x) -1*(w_maxHR_25' * x),x0,A,b,Aeq,beq,lb,ub);
av_max = xhat_max_MHR_25' * w_maxHR_25;
Xmax_MHR_25 = 0.5 * (xhat_max_MHR_25 .* (X_u - X_l) + (X_u + X_l));

xhat_max_MHR_30 = fmincon(@(x) -1*(w_maxHR_30' * x),x0,A,b,Aeq,beq,lb,ub);
av_max = xhat_max_MHR_30' * w_maxHR_30;
Xmax_MHR_30 = 0.5 * (xhat_max_MHR_30 .* (X_u - X_l) + (X_u + X_l));

xhat_max_MHR_35 = fmincon(@(x) -1*(w_maxHR_35' * x),x0,A,b,Aeq,beq,lb,ub);
av_max = xhat_max_MHR_35' * w_maxHR_35;
Xmax_MHR_35 = 0.5 * (xhat_max_MHR_35 .* (X_u - X_l) + (X_u + X_l));

xhat_max_MHR_40 = fmincon(@(x) -1*(w_maxHR_40' * x),x0,A,b,Aeq,beq,lb,ub);
av_max = xhat_max_MHR_40' * w_maxHR_40;
Xmax_MHR_40 = 0.5 * (xhat_max_MHR_40 .* (X_u - X_l) + (X_u + X_l));

% min IHR
xhat_min_IHR_25 = fmincon(@(x) 1*(w_intHR_25' * x),x0,A,b,Aeq,beq,lb,ub);
av_min = xhat_min_IHR_25' * w_intHR_25;
Xmin_IHR_25 = 0.5 * (xhat_min_IHR_25 .* (X_u - X_l) + (X_u + X_l));

xhat_min_IHR_30 = fmincon(@(x) 1*(w_intHR_30' * x),x0,A,b,Aeq,beq,lb,ub);
av_min = xhat_min_IHR_30' * w_intHR_30;
Xmin_IHR_30 = 0.5 * (xhat_min_IHR_30 .* (X_u - X_l) + (X_u + X_l));

xhat_min_IHR_35 = fmincon(@(x) 1*(w_intHR_35' * x),x0,A,b,Aeq,beq,lb,ub);
av_min = xhat_min_IHR_35' * w_intHR_35;
Xmin_IHR_35 = 0.5 * (xhat_min_IHR_35 .* (X_u - X_l) + (X_u + X_l));

xhat_min_IHR_40 = fmincon(@(x) 1*(w_intHR_40' * x),x0,A,b,Aeq,beq,lb,ub);
av_min = xhat_min_IHR_40' * w_intHR_40;
Xmin_IHR_40 = 0.5 * (xhat_min_IHR_40 .* (X_u - X_l) + (X_u + X_l));

% max IHR
xhat_max_IHR_25 = fmincon(@(x) -1*(w_intHR_25' * x),x0,A,b,Aeq,beq,lb,ub);
av_max = xhat_max_IHR_25' * w_intHR_25;
Xmax_IHR_25 = 0.5 * (xhat_max_IHR_25 .* (X_u - X_l) + (X_u + X_l));

xhat_max_IHR_30 = fmincon(@(x) -1*(w_intHR_30' * x),x0,A,b,Aeq,beq,lb,ub);
av_max = xhat_max_IHR_30' * w_intHR_30;
Xmax_IHR_30 = 0.5 * (xhat_max_IHR_30 .* (X_u - X_l) + (X_u + X_l));

xhat_max_IHR_35 = fmincon(@(x) -1*(w_intHR_35' * x),x0,A,b,Aeq,beq,lb,ub);
av_max = xhat_max_IHR_35' * w_intHR_35;
Xmax_IHR_35 = 0.5 * (xhat_max_IHR_35 .* (X_u - X_l) + (X_u + X_l));

xhat_max_IHR_40 = fmincon(@(x) -1*(w_intHR_40' * x),x0,A,b,Aeq,beq,lb,ub);
av_max = xhat_max_IHR_40' * w_intHR_40;
Xmax_IHR_40 = 0.5 * (xhat_max_IHR_40 .* (X_u - X_l) + (X_u + X_l));

%% uncertainty factor sample values (multiply the nominal)
disp('UF sample values:');disp('=================');
disp(' > max(HR) 25');
disp(Xmax_MHR_25);
disp(' > max(int(HR)) 25');
disp(Xmax_IHR_25);
disp(' > max(HR) 30');
disp(Xmax_MHR_30);
disp(' > max(int(HR)) 30');
disp(Xmax_IHR_30);
disp(' > max(HR) 35');
disp(Xmax_MHR_35);
disp(' > max(int(HR)) 35');
disp(Xmax_IHR_35);
%%
FFR = 4;
with_minmax = true;
switch FFR
    case 1
        w_maxHR = w_maxHR_25;
        w_intHR = w_intHR_25;
        xhat_min_MHR = xhat_min_MHR_25;
        xhat_max_MHR = xhat_max_MHR_25;
        xhat_min_IHR = xhat_min_IHR_25;
        xhat_max_IHR = xhat_max_IHR_25;
    case 2
        w_maxHR = w_maxHR_30;
        w_intHR = w_intHR_30;
        xhat_min_MHR = xhat_min_MHR_30;
        xhat_max_MHR = xhat_max_MHR_30;
        xhat_min_IHR = xhat_min_IHR_30;
        xhat_max_IHR = xhat_max_IHR_30;
    case 3
        w_maxHR = w_maxHR_35;
        w_intHR = w_intHR_35;
        xhat_min_MHR = xhat_min_MHR_35;
        xhat_max_MHR = xhat_max_MHR_35;
        xhat_min_IHR = xhat_min_IHR_35;
        xhat_max_IHR = xhat_max_IHR_35;
    case 4
        w_maxHR = w_maxHR_40;
        w_intHR = w_intHR_40;
        xhat_min_MHR = xhat_min_MHR_40;
        xhat_max_MHR = xhat_max_MHR_40;
        xhat_min_IHR = xhat_min_IHR_40;
        xhat_max_IHR = xhat_max_IHR_40;
end        

% sufficient summary plot (maxHR)
figure(1)
plot(Xhat' * w_maxHR, MAX_HR(1:M,FFR), 'o',...
    'MarkerSize', 10,...
    'LineWidth', 1,...
    'MarkerEdgeColor', [0.5, 0.5, 0.5],...
    'MarkerFaceColor', [0.5, 0.5, 0.5]);
if with_minmax
    hold on;
    plot([xhat_min_MHR,xhat_max_MHR]' * w_maxHR, MAX_HR(M+1:M+2,FFR), 'o',...
        'MarkerSize', 10,...
        'LineWidth', 1,...
        'MarkerEdgeColor', [0.9, 0.5, 0.5],...
        'MarkerFaceColor', [0.9, 0.5, 0.5]);
    plot([xhat_min_IHR,xhat_max_IHR]' * w_maxHR, MAX_HR(M+3:M+4,FFR), 'o',...
        'MarkerSize', 10,...
        'LineWidth', 2,...
        'MarkerEdgeColor', [0.9, 0.5, 0.5],...
        'MarkerFaceColor', [1, 1, 1]);
    hold off;
end
ax = gca;
ax.FontSize = 18;
ax.FontName = 'Times New Roman';
ax.Box = 'off';
ax.TickDir = 'out';

% sufficient summary plot (intHR)
figure(2)
plot(Xhat' * w_intHR, INT_HR(1:M,FFR), 'o',...
    'MarkerSize', 10,...
    'LineWidth', 1,...
    'MarkerEdgeColor', [0.5, 0.5, 0.5],...
    'MarkerFaceColor', [0.5, 0.5, 0.5]);
if with_minmax
    hold on;
    plot([xhat_min_MHR,xhat_max_MHR]' * w_intHR, INT_HR(M+1:M+2,FFR), 'o',...
        'MarkerSize', 10,...
        'LineWidth', 1,...
        'MarkerEdgeColor', [0.9, 0.5, 0.5],...
        'MarkerFaceColor', [0.9, 0.5, 0.5]);
    plot([xhat_min_IHR,xhat_max_IHR]' * w_intHR, INT_HR(M+3:M+4,FFR), 'o',...
        'MarkerSize', 10,...
        'LineWidth', 2,...
        'MarkerEdgeColor', [0.9, 0.5, 0.5],...
        'MarkerFaceColor', [1, 1, 1]);
    hold off;
end
ax = gca;
ax.FontSize = 18;
ax.FontName = 'Times New Roman';
ax.Box = 'off';
ax.TickDir = 'out';

%%
% bar chart evaluating active variable
figure(3)
b = bar(cat(2,w_maxHR_25,w_maxHR_30,w_maxHR_35,w_maxHR_40),1,...
            'EdgeColor', 'w');
        
b(1).FaceColor = [0.5, 0.5, 0.5];
b(2).FaceColor = [24,87,155]./255;
b(3).FaceColor = [155,24,87]./255;
b(4).FaceColor = [87,155,24]./255;
ax = gca;
ax.FontSize = 20;
ax.FontName = 'Times New Roman';
ax.Box = 'off';
ax.TickDir = 'out';
ax.XTickLabelRotation = 0;
fig3 = gcf;
fig3.OuterPosition(3) = fig3.OuterPosition(3) + 50;

% bar chart evaluating active variable
figure(4)
b = bar(cat(2,w_intHR_25,w_intHR_30,w_intHR_35,w_intHR_40),1,...
            'EdgeColor', 'w');
        
b(1).FaceColor = [0.5, 0.5, 0.5];
b(2).FaceColor = [24,87,155]./255;
b(3).FaceColor = [155,24,87]./255;
b(4).FaceColor = [87,155,24]./255;
ax = gca;
ax.FontSize = 20;
ax.FontName = 'Times New Roman';
ax.Box = 'off';
ax.TickDir = 'out';
ax.XTickLabelRotation = 0;
fig4 = gcf;
fig4.OuterPosition(3) = fig4.OuterPosition(3) + 50;

%% bar chart evaluating max sample
% MHR cmparison
figure(5)
b = bar(cat(2,xhat_max_MHR_25,xhat_max_MHR_30,xhat_max_MHR_35,xhat_max_MHR_40),1,...
            'EdgeColor', 'w');
b(1).FaceColor = [0.5, 0.5, 0.5];
b(2).FaceColor = [24,87,155]./255;
b(3).FaceColor = [155,24,87]./255;
b(4).FaceColor = [87,155,24]./255;
ax = gca;
ax.FontSize = 20;
ax.FontName = 'Times New Roman';
ax.Box = 'off';
ax.TickDir = 'out';
ax.XTickLabelRotation = 0;
fig5 = gcf;
fig5.OuterPosition(3) = fig5.OuterPosition(3) + 50;

% IHR comparison
figure(6)
b = bar(cat(2,xhat_max_IHR_25,xhat_max_IHR_30,xhat_max_IHR_35,xhat_max_IHR_40),1,...
            'EdgeColor', 'w');
        
b(1).FaceColor = [0.5, 0.5, 0.5];
b(2).FaceColor = [24,87,155]./255;
b(3).FaceColor = [155,24,87]./255;
b(4).FaceColor = [87,155,24]./255;
ax = gca;
ax.FontSize = 20;
ax.FontName = 'Times New Roman';
ax.Box = 'off';
ax.TickDir = 'out';
ax.XTickLabelRotation = 0;
fig6 = gcf;
fig6.OuterPosition(3) = fig6.OuterPosition(3) + 50;

% MHR vs IHR
figure(7)
b = bar(cat(2,xhat_max_MHR_25,xhat_max_IHR_25),1,...
            'EdgeColor', 'w');
        
b(1).FaceColor = [0.5, 0.5, 0.5];
b(2).FaceColor = [24,87,155]./255;
ax = gca;
ax.FontSize = 20;
ax.FontName = 'Times New Roman';
ax.Box = 'off';
ax.TickDir = 'out';
ax.XTickLabelRotation = 0;
fig7 = gcf;
fig7.OuterPosition(3) = fig7.OuterPosition(3) + 50;
