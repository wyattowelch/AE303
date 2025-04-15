%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AE302 Lab 3 - Wyatt Welch
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear all, load('Raw4.mat'), load('Raw5.mat'), load('Raw6.mat'), load('TP4.mat'), load('TP5.mat'), load('TP6.mat')

% Data Setup
Raw4 = mean(table2array(Raw4),2);
Raw5 = mean(table2array(Raw5),2);
Raw6 = mean(table2array(Raw6),2);

for i = 1:size(Raw4, 1)-4
    Raw4(i+4, :) = Raw4(i+4, :) - Raw4(mod(i-1,4)+1, :); 
end
Rows = [1:4,41:44];
rows = [1,11];
Raw4 = Raw4(~ismember(1:size(Raw4, 1), Rows), :);
Raw4 = mean(Raw4, 2); 
Raw4 = Raw4 * 6894.76;

for i = 1:size(Raw5, 1)-4
    Raw5(i+4, :) = Raw5(i+4, :) - Raw5(mod(i-1,4)+1, :); 
end
Raw5 = Raw5(~ismember(1:size(Raw5, 1), Rows), :);
Raw5 = mean(Raw5, 2); 
Raw5 = Raw5 * 6894.76;

for i = 1:size(Raw6, 1)-4
    Raw6(i+4, :) = Raw6(i+4, :) - Raw6(mod(i-1,4)+1, :); 
end
Raw6 = Raw6(~ismember(1:size(Raw6, 1), Rows), :);
Raw6 = mean(Raw6, 2); 
Raw6 = Raw6 * 6894.76;

TP4 = table2array(TP4);
TP4 = TP4(~ismember(1:size(TP4, 1), rows), :);
TP5 = table2array(TP5);
TP5 = TP5(~ismember(1:size(TP5, 1), rows), :);
TP6 = table2array(TP6);
TP6 = TP6(~ismember(1:size(TP6, 1), rows), :);

T_amb4 = (76.6 - 32) * (5/9) + 273.15;
T_amb5 = (73.6 - 32) * (5/9) + 273.15;
T_amb6 = (75.7 - 32) * (5/9) + 273.15;

T4 = (TP4(:,2) - 32) .* (5/9) + 273.15;
T5 = (TP5(:,2) - 32) .* (5/9) + 273.15;
T6 = (TP6(:,2) - 32) .* (5/9) + 273.15;

P_a4 = 29.94 * 3386.39;
P_a5 = 29.93 * 3386.39;
P_a6 = 29.93 * 3386.39;

Raw4 = Raw4 + P_a4;
Raw5 = Raw5 + P_a5;
Raw6 = Raw6 + P_a6;

P4 = TP4(:,1) .* 249.0889;
P5 = TP5(:,1) .* 249.0889;
P6 = TP6(:,1) .* 249.0889;

R = 287.05;

D4 = 0.1016;
D5 = 0.1267;
D6 = 0.1524;

% Calculations
% 1 = Static, 2 = Total, 3 = Rear, 4 = Forward

delP4 = Raw4(4:4:end,:) - Raw4(3:4:end,:);
delP5 = Raw5(4:4:end,:) - Raw5(3:4:end,:);
delP6 = Raw6(4:4:end,:) - Raw6(3:4:end,:);

q4 = Raw4(2:4:end,:) - Raw4(1:4:end,:);
q5 = Raw5(2:4:end,:) - Raw5(1:4:end,:);
q6 = Raw6(2:4:end,:) - Raw6(1:4:end,:);

dP_div_q4 = delP4 ./ q4;
dP_div_q5 = delP5 ./ q5;
dP_div_q6 = delP6 ./ q6;
dP_div_q_crit = 1.220;

den4 = abs(Raw4(1:4:end,:)) ./ (R .* T4);
den5 = abs(Raw5(1:4:end,:)) ./ (R .* T5);
den6 = abs(Raw6(1:4:end,:)) ./ (R .* T6);

Visc4 = (1.716e-5) .* ((T4 ./ 273.15) .^ 1.5) .* ((273.15 + 110.4) ./ (T4 + 110.4));
Visc5 = (1.716e-5) .* ((T5 ./ 273.15) .^ 1.5) .* ((273.15 + 110.4) ./ (T5 + 110.4));
Visc6 = (1.716e-5) .* ((T6 ./ 273.15) .^ 1.5) .* ((273.15 + 110.4) ./ (T6 + 110.4));

U4 = sqrt((2 .* q4) ./ den4);
U5 = sqrt((2 .* q5) ./ den5);
U6 = sqrt((2 .* q6) ./ den6);

Re4 = (den4 .* U4 .* D4) ./ Visc4;
Re5 = (den5 .* U5 .* D5) ./ Visc5;
Re6 = (den6 .* U6 .* D6) ./ Visc6;

Re4_interp = interp1(dP_div_q4, Re4, 1.220); 
Re5_interp = interp1(dP_div_q5, Re5, 1.220); 
Re6_interp = interp1(dP_div_q6, Re6, 1.220); 
Re_crit = [Re4_interp Re5_interp Re6_interp];

TF4 = 3.85e5 / Re4_interp;
TF5 = 3.85e5 / Re5_interp;
TF6 = 3.85e5 / Re6_interp;
TF = [TF4 TF5 TF6];

Tu4 = (TF4 - 1) / 0.15;
Tu5 = (TF5 - 1) / 0.15;
Tu6 = (TF6 - 1) / 0.15;
Tu = [Tu4 Tu5 Tu6];

Re_eff4 = TF4 * Re4;
Re_eff5 = TF5 * Re5;
Re_eff6 = TF6 * Re6;


% Plots

figure(1)
hold on, grid on
plot(Re4, dP_div_q4, '-o', 'Color', 'red')
plot(Re5, dP_div_q5, '-o', 'Color', 'blue')
plot(Re6, dP_div_q6, '-o', 'Color', 'magenta')
plot(Re4_interp, 1.220, 'kx', 'MarkerSize', 10, 'LineWidth', 2)
plot(Re5_interp, 1.220, 'kx', 'MarkerSize', 10, 'LineWidth', 2)
plot(Re6_interp, 1.220, 'kx', 'MarkerSize', 10, 'LineWidth', 2)
plot([min([Re4; Re5; Re6]), max([Re4; Re5; Re6])], [1.220, 1.220], 'Color', 'black', 'LineStyle', '--', 'LineWidth', 0.5);

xlabel('Reynolds Number (Re)', 'Interpreter', 'latex')
ylabel('$\frac{\Delta P}{q}$', 'Interpreter', 'latex')
title('$\Delta P / q$ vs Reynolds Number', 'Interpreter', 'latex')
criticalReString = sprintf('Critical Re (%.2e, %.2e, %.2e)', Re4_interp, Re5_interp, Re6_interp);
legend('4" Sphere', '4.987" Sphere', '6" Sphere', criticalReString, 'Location', 'best')


figure(2)
hold on, grid on
yyaxis left;
plot(Re4_interp, TF4, 'o', 'Color', 'blue', 'DisplayName', '4" Sphere (TF)');
plot(Re5_interp, TF5, 'o', 'Color', 'green', 'DisplayName', '4.987" Sphere (TF)');
plot(Re6_interp, TF6, 'o', 'Color', 'magenta', 'DisplayName', '6" Sphere (TF)');
ylabel('Turbulence Factor (TF)');

yyaxis right;
plot(Re4_interp, Tu4, 's', 'Color', 'blue', 'DisplayName', '4" Sphere (Tu)');
plot(Re5_interp, Tu5, 's', 'Color', 'green', 'DisplayName', '4.987" Sphere (Tu)');
plot(Re6_interp, Tu6, 's', 'Color', 'magenta', 'DisplayName', '6" Sphere (Tu)');
ylabel('Percent Turbulence (Tu)');

xlabel('Critical Unit Reynolds Number (Re_{critical})');
legend('Location', 'northeast');



fprintf("Completed Run \n")
