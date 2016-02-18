% active subspace post-processing
clear all; clc;

M = 24;
m = 12;

MAX_HR = zeros(24,1);
INT_HR = zeros(24,1);
MAX_HR_N = zeros(24,1);
INT_HR_N = zeros(24,1);

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
];

MAX_HR_N(:,1) = [
    1.82043E+11
    1.78185E+11
    2.03071E+11
    2.16750E+11
    1.68494E+11
    2.07378E+11
    1.67374E+11
    1.62586E+11
    2.21614E+11
    2.11214E+11
    1.73034E+11
    2.21004E+11
    1.96252E+11
    2.03434E+11
    2.13360E+11
    2.14432E+11
    2.11036E+11
    2.15261E+11
    1.85174E+11
    1.94482E+11
    1.61116E+11
    2.10547E+11
    2.24437E+11
    2.12250E+11
];

INT_HR_N(:,1) = [
    2.165713E+06
    2.184937E+06
    2.354858E+06
    2.342367E+06
    2.103353E+06
    2.236647E+06
    2.047337E+06
    1.981653E+06
    2.436046E+06
    2.332522E+06
    2.027786E+06
    2.417745E+06
    2.274784E+06
    2.242562E+06
    2.435799E+06
    2.407126E+06
    2.420863E+06
    2.384889E+06
    2.263310E+06
    2.251954E+06
    1.916489E+06
    2.257475E+06
    2.539546E+06
    2.323751E+06
];


%% Normalized Data
% Xhat: samples normalized by min/max values [-1,1]

% X: raw sample values
% X_u: vector of max sample values (upper)
% X_l: vector of min sample values (lower)
load /Users/memory/Dropbox/Matlab/Hyshot_QMU/chem24_samples.mat

Xhat = (2 * X - repmat(X_u + X_l,1,M)) ./ repmat(X_u - X_l,1,M);

%% compute linear regression weights
close all;

uhat = [ones(M,1) Xhat'] \ MAX_HR;
w_maxHR = uhat(2:m+1) / norm(uhat(2:m+1));

uhat = [ones(M,1) Xhat'] \ INT_HR;
w_intHR = uhat(2:m+1) / norm(uhat(2:m+1));

disp('weights (maxHR):')
disp(w_maxHR)

disp('weights (intHR):')
disp(w_intHR)
    
A = [];
b = [];
Aeq = [];
beq = [];
lb = -1 * ones(m,1);
ub = 1* ones(m,1);
x0 = zeros(m,1);

x_min = fmincon(@(x) w_maxHR' * x,x0,A,b,Aeq,beq,lb,ub);
av_min = x_min' * w_maxHR;
disp('maxHR min a.v.:')
disp(av_min)
disp('maxHR min weights:')
disp(x_min')

x_max = fmincon(@(x) -1*(w_maxHR' * x),x0,A,b,Aeq,beq,lb,ub);
av_max = x_max' * w_maxHR;
disp('maxHR max a.v.:')
disp(av_max)
disp('maxHR max weights:')
disp(x_max')
Xmax_MHR = 0.5 * (x_max .* (X_u - X_l) + (X_u + X_l));

%%

x_min = fmincon(@(x) w_intHR' * x,x0,A,b,Aeq,beq,lb,ub);
av_min = x_min' * w_intHR;
disp('intHR min a.v.:')
disp(av_min)
disp('intHR min weights:')
disp(x_min')

x_max = fmincon(@(x) -1*(w_intHR' * x),x0,A,b,Aeq,beq,lb,ub);
av_max = x_max' * w_intHR;
disp('intHR max a.v.:')
disp(av_max)
disp('intHR max weights:')
disp(x_max')
Xmax_IHR = 0.5 * (x_max .* (X_u - X_l) + (X_u + X_l));

%%
% for plotting
FFR = 1;

% sufficient summary plot (maxHR)
figure(1)
plot(Xhat' * w_maxHR, MAX_HR, 'o',...
    'MarkerSize', 10,...
    'LineWidth', 1,...
    'MarkerEdgeColor', [0.5, 0.5, 0.5],...
    'MarkerFaceColor', [0.5, 0.5, 0.5]);
ax = gca;
ax.FontSize = 18;
ax.FontName = 'Times New Roman';
ax.Box = 'off';
ax.TickDir = 'out';

% sufficient summary plot (intHR)
figure(2)
plot(Xhat' * w_intHR, INT_HR, 'o',...
    'MarkerSize', 10,...
    'LineWidth', 1,...
    'MarkerEdgeColor', [0.5, 0.5, 0.5],...
    'MarkerFaceColor', [0.5, 0.5, 0.5]);
ax = gca;
ax.FontSize = 18;
ax.FontName = 'Times New Roman';
ax.Box = 'off';
ax.TickDir = 'out';


% bar chart evaluating weights
figure(3)
b = bar(cat(2,w_maxHR,w_intHR),1,...
            'EdgeColor', 'w');
        
b(1).FaceColor = [0.5, 0.5, 0.5];
b(2).FaceColor = [24,87,155]./255;
ax = gca;
ax.FontSize = 20;
ax.FontName = 'Times New Roman';
ax.Box = 'off';
ax.TickDir = 'out';
ax.XTickLabelRotation = 0;
fig3 = gcf;
fig3.OuterPosition(3) = fig3.OuterPosition(3) + 50;





