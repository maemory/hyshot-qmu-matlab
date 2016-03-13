% active subspace post-processing
clear all; clc;

M = 24;
m = 12;

MAX_HR = zeros(M+4,1);
INT_HR = zeros(M+4,1);
EP = zeros(M+4,1);
SV = zeros(M+4,1);

MAX_HR(:,1) = [
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

INT_HR(:,1) = [
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

EP(:,1) = [
    2.389668
    2.393329
    2.402416
    2.399836
    2.386293
    2.388704
    2.389445
    2.385821
    2.403574
    2.393145
    2.380624
    2.400453
    2.394596
    2.388582
    2.408596
    2.405758
    2.401817
    2.398059
    2.397976
    2.400359
    2.379926
    2.393892
    2.410452
    2.400368
    2.369446
    2.404313
    0
    2.404178
];

SV(:,1) = [
    98.84409
    98.82307
    98.80643
    98.80103
    98.84968
    98.83205
    98.84048
    98.85662
    98.77254
    98.81339
    98.86236
    98.79408
    98.82126
    98.82803
    98.77375
    98.78320
    98.78914
    98.79341
    98.81894
    98.80613
    98.86254
    98.82293
    98.75828
    98.80117
    99.00358
    98.74473
    0
    98.75498
];
SV = 100-SV;


%% Normalized Data
% Xhat: samples normalized by min/max values [-1,1]

% X: raw sample values
% X_u: vector of max sample values (upper)
% X_l: vector of min sample values (lower)
load /Users/memory/Dropbox/Matlab/Hyshot_QMU/PAUL-chem24_samples.mat

Xhat = (2 * X - repmat(X_u + X_l,1,M)) ./ repmat(X_u - X_l,1,M);

%% compute linear regression weights
close all;

uhat = [ones(M,1) Xhat'] \ MAX_HR(1:M);
w_maxHR = uhat(2:m+1) / norm(uhat(2:m+1));

uhat = [ones(M,1) Xhat'] \ INT_HR(1:M);
w_intHR = uhat(2:m+1) / norm(uhat(2:m+1));

uhat = [ones(M,1) Xhat'] \ EP(1:M);
w_EP = uhat(2:m+1) / norm(uhat(2:m+1));

uhat = [ones(M,1) Xhat'] \ SV(1:M);
w_SV = uhat(2:m+1) / norm(uhat(2:m+1));

disp('weights (maxHR):')
disp(w_maxHR)

disp('weights (intHR):')
disp(w_intHR)

disp('weights (exit press):')
disp(w_EP)

disp('weights (subsonic vol):')
disp(w_SV)
    
A = [];
b = [];
Aeq = [];
beq = [];
lb = -1 * ones(m,1);
ub = 1* ones(m,1);
x0 = zeros(m,1);

%% maximum HR
xhat_min_MHR = fmincon(@(x) w_maxHR' * x,x0,A,b,Aeq,beq,lb,ub);
av_min = xhat_min_MHR' * w_maxHR;
disp('maxHR min a.v.:')
disp(av_min)
disp('maxHR min weights:')
disp(xhat_min_MHR')
Xmin_MHR = 0.5 * (xhat_min_MHR .* (X_u - X_l) + (X_u + X_l));

xhat_max_MHR = fmincon(@(x) -1*(w_maxHR' * x),x0,A,b,Aeq,beq,lb,ub);
av_max = xhat_max_MHR' * w_maxHR;
disp('maxHR max a.v.:')
disp(av_max)
disp('maxHR max weights:')
disp(xhat_max_MHR')
Xmax_MHR = 0.5 * (xhat_max_MHR .* (X_u - X_l) + (X_u + X_l));

%% integrated HR

xhat_min_IHR = fmincon(@(x) w_intHR' * x,x0,A,b,Aeq,beq,lb,ub);
av_min = xhat_min_IHR' * w_intHR;
disp('intHR min a.v.:')
disp(av_min)
disp('intHR min weights:')
disp(xhat_min_IHR')
Xmin_IHR = 0.5 * (xhat_min_IHR .* (X_u - X_l) + (X_u + X_l));


xhat_max_IHR = fmincon(@(x) -1*(w_intHR' * x),x0,A,b,Aeq,beq,lb,ub);
av_max = xhat_max_IHR' * w_intHR;
disp('intHR max a.v.:')
disp(av_max)
disp('intHR max weights:')
disp(xhat_max_IHR')
Xmax_IHR = 0.5 * (xhat_max_IHR .* (X_u - X_l) + (X_u + X_l));

%% exit pressure

xhat_min_EP = fmincon(@(x) w_EP' * x,x0,A,b,Aeq,beq,lb,ub);
av_min = xhat_min_EP' * w_EP;
disp('EP min a.v.:')
disp(av_min)
disp('EP min weights:')
disp(xhat_min_EP')
Xmin_EP = 0.5 * (xhat_min_EP .* (X_u - X_l) + (X_u + X_l));


xhat_max_EP = fmincon(@(x) -1*(w_EP' * x),x0,A,b,Aeq,beq,lb,ub);
av_max = xhat_max_EP' * w_EP;
disp('EP max a.v.:')
disp(av_max)
disp('EP max weights:')
disp(xhat_max_EP')
Xmax_EP = 0.5 * (xhat_max_EP .* (X_u - X_l) + (X_u + X_l));

%% subsonic volume

xhat_min_SV = fmincon(@(x) w_SV' * x,x0,A,b,Aeq,beq,lb,ub);
av_min = xhat_min_SV' * w_SV;
disp('SV min a.v.:')
disp(av_min)
disp('SV min weights:')
disp(xhat_min_SV')
Xmin_SV = 0.5 * (xhat_min_SV .* (X_u - X_l) + (X_u + X_l));


xhat_max_SV = fmincon(@(x) -1*(w_SV' * x),x0,A,b,Aeq,beq,lb,ub);
av_max = xhat_max_SV' * w_SV;
disp('SV max a.v.:')
disp(av_max)
disp('SV max weights:')
disp(xhat_max_SV')
Xmax_SV = 0.5 * (xhat_max_SV .* (X_u - X_l) + (X_u + X_l));

%% 

X_labels={
 'r1',
 'r2',
 'r5',
 'r8',
 'r9/r10',
 'r11',
 'r14/r15',
 'r16',
 'r17',
 'r19',
 'r22',
 'r24'
};

for i=1:12
   disp(sprintf('%f %f %s',X_l(i),X_u(i),X_labels{i})); 
end

% script-to-Hong reaction
A_mapping = {
'r1',1,
'r2',2,
'r3',2,
'r4',2,
'r5',3,
'r6',4,
'r7',4,
'r8',5,
'r9',6,
'r10',6,
'r11',7,
'r12',7,
'r13',8,
'r14',9,
'r15',9,
'r16',10,
'r17',11,
'r18',12,
'r19',13,
'r20',14,
'r21',15,
'r22',16,
'r23',17,
'r24',18,
'r25',18,
'r26',18,
'r27',18,
'r28',19,
'r29',20
};

%% uncertainty factor sample values (multiply the nominal)
disp('UF sample values:');disp('=================');
disp(' > min(HR)');
disp(Xmin_MHR);
disp(' > max(HR)');
disp(Xmax_MHR);
disp(' > min(int(HR))');
disp(Xmin_IHR);
disp(' > max(int(HR))');
disp(Xmax_IHR);

%%

% sufficient summary plot (maxHR)
figure(1)
plot(Xhat' * w_maxHR, MAX_HR(1:M), 'o',...
    'MarkerSize', 10,...
    'LineWidth', 1,...
    'MarkerEdgeColor', [0.5, 0.5, 0.5],...
    'MarkerFaceColor', [0.5, 0.5, 0.5]);
hold on;
plot([xhat_min_MHR,xhat_max_MHR]' * w_maxHR, MAX_HR(25:26), 'o',...
        'MarkerSize', 10,...
        'LineWidth', 1,...
        'MarkerEdgeColor', [0.9, 0.5, 0.5],...
        'MarkerFaceColor', [0.9, 0.5, 0.5]);
% plot([xhat_min_IHR,xhat_max_IHR]' * w_maxHR, MAX_HR(27:28), 'o',...
%         'MarkerSize', 10,...
%         'LineWidth', 2,...
%         'MarkerEdgeColor', [0.9, 0.5, 0.5],...
%         'MarkerFaceColor', [1, 1, 1]);
hold off;
ax = gca;
ax.FontSize = 18;
ax.FontName = 'Times New Roman';
ax.Box = 'off';
ax.TickDir = 'out';

% sufficient summary plot (intHR)
figure(2)
plot(Xhat' * w_intHR, INT_HR(1:M), 'o',...
    'MarkerSize', 10,...
    'LineWidth', 1,...
    'MarkerEdgeColor', [0.5, 0.5, 0.5],...
    'MarkerFaceColor', [0.5, 0.5, 0.5]);
hold on;
plot([xhat_min_MHR,xhat_max_MHR]' * w_intHR, INT_HR(25:26), 'o',...
        'MarkerSize', 10,...
        'LineWidth', 1,...
        'MarkerEdgeColor', [0.9, 0.5, 0.5],...
        'MarkerFaceColor', [0.9, 0.5, 0.5]);
% plot([xhat_min_IHR,xhat_max_IHR]' * w_intHR, INT_HR(27:28), 'o',...
%         'MarkerSize', 10,...
%         'LineWidth', 2,...
%         'MarkerEdgeColor', [0.9, 0.5, 0.5],...
%         'MarkerFaceColor', [1, 1, 1]);
hold off;
ax = gca;
ax.FontSize = 18;
ax.FontName = 'Times New Roman';
ax.Box = 'off';
ax.TickDir = 'out';

%%
% sufficient summary plot (EP)
figure(1)
plot(Xhat' * w_EP, EP(1:M), 'o',...
    'MarkerSize', 10,...
    'LineWidth', 1,...
    'MarkerEdgeColor', [0.5, 0.5, 0.5],...
    'MarkerFaceColor', [0.5, 0.5, 0.5]);
hold on;
plot([xhat_min_MHR,xhat_max_MHR]' * w_EP, EP(25:26), 'o',...
        'MarkerSize', 10,...
        'LineWidth', 1,...
        'MarkerEdgeColor', [0.9, 0.5, 0.5],...
        'MarkerFaceColor', [0.9, 0.5, 0.5]);
plot([xhat_max_IHR]' * w_EP, EP(28), 'o',...
        'MarkerSize', 10,...
        'LineWidth', 2,...
        'MarkerEdgeColor', [0.9, 0.5, 0.5],...
        'MarkerFaceColor', [1, 1, 1]);
hold off;
ax = gca;
ax.FontSize = 18;
ax.FontName = 'Times New Roman';
ax.Box = 'off';
ax.TickDir = 'out';

% sufficient summary plot (SV)
figure(2)
plot(Xhat' * w_SV, SV(1:M), 'o',...
    'MarkerSize', 10,...
    'LineWidth', 1,...
    'MarkerEdgeColor', [0.5, 0.5, 0.5],...
    'MarkerFaceColor', [0.5, 0.5, 0.5]);
hold on;
plot([xhat_min_MHR,xhat_max_MHR]' * w_SV, SV(25:26), 'o',...
        'MarkerSize', 10,...
        'LineWidth', 1,...
        'MarkerEdgeColor', [0.9, 0.5, 0.5],...
        'MarkerFaceColor', [0.9, 0.5, 0.5]);
plot([xhat_max_IHR]' * w_SV, SV(28), 'o',...
        'MarkerSize', 10,...
        'LineWidth', 2,...
        'MarkerEdgeColor', [0.9, 0.5, 0.5],...
        'MarkerFaceColor', [1, 1, 1]);
hold off;
ax = gca;
ax.FontSize = 18;
ax.FontName = 'Times New Roman';
ax.Box = 'off';
ax.TickDir = 'out';

%% bar chart evaluating active variable
figure(3)
b = bar(cat(2,w_maxHR,w_intHR,w_EP,w_SV),1,...
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

%% bar chart evaluating weights
figure(4)
b = bar(cat(2,xhat_min_MHR,xhat_min_IHR,xhat_min_EP,xhat_min_SV),1,...
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

figure(5)
b = bar(cat(2,xhat_max_MHR,xhat_max_IHR,xhat_max_EP,xhat_max_SV),1,...
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

