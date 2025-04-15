%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AE302 Lab 4 - Wyatt Welch
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear all

% Load Data & Convert to Tables in Pa

load("q0a0.mat"); load("q5aN5.mat"); load("q5a0.mat"); load("q5a5.mat"); 
load("q5a10.mat"); load("q5a15.mat"); load("q5a20_1.mat"); load("q5a20_2.mat");

q0a0 = mean(table2array(q0a0),2) .* 6894.76;                q5aN5 = mean(table2array(q5aN5),2) .* 6894.76 - q0a0;    
q5a0 = mean(table2array(q5a0),2) .* 6894.76 - q0a0;         q5a5 = mean(table2array(q5a5),2) .* 6894.76 - q0a0; 
q5a10 = mean(table2array(q5a10),2) .* 6894.76 - q0a0;       q5a15 = mean(table2array(q5a15),2) .* 6894.76 - q0a0; 
q5a20_1 = mean(table2array(q5a20_1),2) .* 6894.76 - q0a0;   q5a20_2 = mean(table2array(q5a20_2),2) .* 6894.76 - q0a0; % Pa
q0a0 = q0a0 - q0a0;

load("Theo_Vals.mat"); Theo = table2array(Theo_Vals);
N5x = Theo(:,1);    N5y = Theo(:,2);    N5Cp = Theo(:,3); 
P0x = Theo(:,4);    P0y = Theo(:,5);    P0Cp = Theo(:,6); 
P5x = Theo(:,7);    P5y = Theo(:,8);    P5Cp = Theo(:,9); 
P10x = Theo(:,10);  P10y = Theo(:,11);  P10Cp = Theo(:,12); 
P15x = Theo(:,13);  P15y = Theo(:,14);  P15Cp = Theo(:,15); 
P20x = Theo(:,16);  P20y = Theo(:,17);  P20Cp = Theo(:,18); 


Pa = 30.11 * 3386.39;               % Pa
Ta = (79.5 - 32) * (5/9) + 273.15;  % K
T_ref = 291.15; 
mu_ref = 1.827e-5;
S = 120;
mu = mu_ref * (Ta/T_ref)^(3/2) * (T_ref + S)/(Ta + S); 

alpha = [-5 0 5 10 15 20];

Xc_lower = [0.015, 0.029, 0.055, 0.080, 0.105, 0.157, 0.207, 0.257, 0.306, ...
    0.407, 0.507, 0.608, 0.708, 0.812, 0.912, 1.000];
Xc_upper = [0.000, 0.013, 0.025, 0.048, 0.073, 0.097, 0.150, 0.200, 0.250, ...
    0.300, 0.400, 0.500, 0.600, 0.700, 0.800, 0.900];
Yc_lower = [-0.0089, -0.0156, -0.0160, -0.0190, -0.0216, -0.0262, -0.0299, ...
    -0.0325, -0.0353, -0.0393, -0.0405, -0.0372, -0.0355, -0.0312, -0.0136, 0.0000];
Yc_upper = [0.0000, 0.0394, 0.0514, 0.0678, 0.0794, 0.0868, 0.0933, 0.0927, ...
    0.0922, 0.0857, 0.0777, 0.0656, 0.0522, 0.0372, 0.0244, 0.0125];
xac_c = 0.238;
yac_c = 0.07;
c = 1; 

Re = 827000;
R = 287;
rho = Pa / (R * Ta);

% Coefficent Calculations

qinf_5N5 = q5aN5(2) - q5aN5(1);  qinf_50 = q5a0(2) - q5a0(1);
qinf_55 = q5a5(2) - q5a5(1);     qinf_510 = q5a10(2) - q5a10(1);
qinf_515 = q5a15(2) - q5a15(1);  qinf_520 = q5a20_1(2) - q5a20_1(1);

taps_5N5 = q5aN5(3:34);  taps_50 = q5a0(3:34); 
taps_55 = q5a5(3:34);    taps_510 = q5a10(3:34); 
taps_515 = q5a15(3:34);  taps_520 = q5a20_1(3:34); 

u_5N5 = sqrt((2 .* (qinf_5N5)) ./ rho);  u_50 = sqrt((2 .* (qinf_50)) ./ rho);
u_55 = sqrt((2 .* (qinf_55)) ./ rho);    u_510 = sqrt((2 .* (qinf_510)) ./ rho);
u_515 = sqrt((2 .* (qinf_515)) ./ rho);  u_520 = sqrt((2 .* (qinf_520)) ./ rho);

Re_5N5 = (rho * u_5N5 * c) / mu;    Re_50 = (rho * u_50 * c) / mu;
Re_55 = (rho * u_55 * c) / mu;       Re_510 = (rho * u_510 * c) / mu;
Re_515 = (rho * u_515 * c) / mu;    Re_520 = (rho * u_520 * c) / mu;
Re = [Re_5N5;Re_50;Re_55;Re_510;Re_515;Re_520];

Cp_5N5 = (taps_5N5 - q5aN5(1)) / qinf_5N5;  Cp_50 = (taps_50 - q5a0(1)) / qinf_50; 
Cp_55 = (taps_55 - q5a5(1)) / qinf_55;      Cp_510 = (taps_510 - q5a10(1)) / qinf_510; 
Cp_515 = (taps_515 - q5a15(1)) / qinf_515;  Cp_520 = (taps_520 - q5a20_1(1)) / qinf_520;

Cp_lower_5N5 = Cp_5N5(1:16);  Cp_lower_50 = Cp_50(1:16); 
Cp_lower_55 = Cp_55(1:16);    Cp_lower_510 = Cp_510(1:16); 
Cp_lower_515 = Cp_515(1:16);  Cp_lower_520 = Cp_520(1:16); 

Cp_upper_5N5 = Cp_5N5(17:32); Cp_upper_50 = Cp_50(17:32); 
Cp_upper_55 = Cp_55(17:32);   Cp_upper_510 = Cp_510(17:32); 
Cp_upper_515 = Cp_515(17:32); Cp_upper_520 = Cp_520(17:32); 

Cn_5N5 = mean(trapz(Xc_lower, Cp_lower_5N5)) - mean(trapz(Xc_upper, Cp_upper_5N5));  
Cn_50 = mean(trapz(Xc_lower, Cp_lower_50)) - mean(trapz(Xc_upper, Cp_upper_50)); 
Cn_55 = mean(trapz(Xc_lower, Cp_lower_55)) - mean(trapz(Xc_upper, Cp_upper_55));   
Cn_510 = mean(trapz(Xc_lower, Cp_lower_510)) - mean(trapz(Xc_upper, Cp_upper_510)); 
Cn_515 = mean(trapz(Xc_lower, Cp_lower_515)) - mean(trapz(Xc_upper, Cp_upper_515));  
Cn_520 = mean(trapz(Xc_lower, Cp_lower_520)) - mean(trapz(Xc_upper, Cp_upper_520)); 

dyu_dx = gradient(Yc_upper) ./ gradient(Xc_upper);
dyl_dx = gradient(Yc_lower) ./ gradient(Xc_lower);

Ca_5N5 = mean(trapz(Xc_upper, Cp_upper_5N5 .* dyu_dx)) - mean(trapz(Xc_lower, Cp_lower_5N5 .* dyl_dx));  
Ca_50 = mean(trapz(Xc_upper, Cp_upper_50 .* dyu_dx)) - mean(trapz(Xc_lower, Cp_lower_50 .* dyl_dx)); 
Ca_55 = mean(trapz(Xc_upper, Cp_upper_55 .* dyu_dx)) - mean(trapz(Xc_lower, Cp_lower_55 .* dyl_dx));   
Ca_510 = mean(trapz(Xc_upper, Cp_upper_510 .* dyu_dx)) - mean(trapz(Xc_lower, Cp_lower_510 .* dyl_dx)); 
Ca_515 = mean(trapz(Xc_upper, Cp_upper_515 .* dyu_dx)) - mean(trapz(Xc_lower, Cp_lower_515 .* dyl_dx));  
Ca_520 = mean(trapz(Xc_upper, Cp_upper_520 .* dyu_dx)) - mean(trapz(Xc_lower, Cp_lower_520 .* dyl_dx)); 

Cl_5N5 = Cn_5N5 * cosd(-5) - Ca_5N5 * sind(-5); Cl_50 = Cn_50 * cosd(0) - Ca_50 * sind(0); 
Cl_55 = Cn_55 * cosd(5) - Ca_55 * sind(5);      Cl_510 = Cn_510 * cosd(10) - Ca_510 * sind(10); 
Cl_515 = Cn_515 * cosd(15) - Ca_515 * sind(15); Cl_520 = Cn_520 * cosd(20) - Ca_520 * sind(20); 
Cl = [Cl_5N5 Cl_50 Cl_55 Cl_510 Cl_515 Cl_520];

Cd_5N5 = Cn_5N5 * sind(-5) + Ca_5N5 * cosd(-5); Cd_50 = Cn_50 * sind(0) + Ca_50 * cosd(0); 
Cd_55 = Cn_55 * sind(5) + Ca_55 * cosd(5);      Cd_510 = Cn_510 * sind(10) + Ca_510 * cosd(10); 
Cd_515 = Cn_515 * sind(15) + Ca_515 * cosd(15); Cd_520 = Cn_520 * sind(20) + Ca_520 * cosd(20); 
Cd = [Cd_5N5 Cd_50 Cd_55 Cd_510 Cd_515 Cd_520];

dyu_dx = dyu_dx';
dyl_dx = dyl_dx';
x_rel = (Xc_upper - xac_c)';
y_rel_u = (Yc_upper - yac_c)';
y_rel_l = (Yc_lower - yac_c)';
Xc_upper = Xc_upper(:); 

t1_5N5 = (Cp_upper_5N5 - Cp_lower_5N5) .* x_rel;    t2_5N5 = Cp_upper_5N5 .* dyu_dx .* y_rel_u;    t3_5N5 = Cp_lower_5N5 .* dyl_dx .* y_rel_l;
t1_50 = (Cp_upper_50 - Cp_lower_50) .* x_rel;       t2_50 = Cp_upper_50 .* dyu_dx .* y_rel_u;      t3_50 = Cp_lower_50 .* dyl_dx .* y_rel_l;
t1_55 = (Cp_upper_55 - Cp_lower_55) .* x_rel;       t2_55 = Cp_upper_55 .* dyu_dx .* y_rel_u;      t3_55 = Cp_lower_55 .* dyl_dx .* y_rel_l;
t1_510 = (Cp_upper_510 - Cp_lower_510) .* x_rel;    t2_510 = Cp_upper_510 .* dyu_dx .* y_rel_u;    t3_510 = Cp_lower_510 .* dyl_dx .* y_rel_l;
t1_515 = (Cp_upper_515 - Cp_lower_515) .* x_rel;    t2_515 = Cp_upper_515 .* dyu_dx .* y_rel_u;    t3_515 = Cp_lower_515 .* dyl_dx .* y_rel_l;
t1_520 = (Cp_upper_520 - Cp_lower_520) .* x_rel;    t2_520 = Cp_upper_520 .* dyu_dx .* y_rel_u;    t3_520 = Cp_lower_520 .* dyl_dx .* y_rel_l;

Cm_ac_5N5 = trapz(Xc_upper, t1_5N5 - t2_5N5 + t3_5N5);  Cm_ac_50 = trapz(Xc_upper, t1_50 - t2_50 + t3_50);
Cm_ac_55 = trapz(Xc_upper, t1_55 - t2_55 + t3_55);      Cm_ac_510 = trapz(Xc_upper, t1_510 - t2_510 + t3_510);
Cm_ac_515 = trapz(Xc_upper, t1_515 - t2_515 + t3_515);  Cm_ac_520 = trapz(Xc_upper, t1_520 - t2_520 + t3_520);
Cm_ac = mean([Cm_ac_5N5; Cm_ac_50; Cm_ac_55; Cm_ac_510; Cm_ac_515; Cm_ac_520],2)';

% Wake Work
yw_in = [0, 2, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10, 10.5, 11, 13];
yw = (yw_in - mean(yw_in)) * .0254;

qw_5N5 = q5aN5(35:54);  qw_50 = q5a0(35:54);
qw_55 = q5a5(35:54);    qw_510 = q5a10(35:54);
qw_515 = q5a15(35:54);  qw_520 = q5a20_1(35:54);

uw_5N5 = sqrt((2 * qw_5N5) / rho);  uw_50 = sqrt((2 * qw_50) / rho);
uw_55 = sqrt((2 * qw_55) / rho);    uw_510 = sqrt((2 * qw_510) / rho);
uw_515 = sqrt((2 * qw_515) / rho);  uw_520 = sqrt((2 * qw_520) / rho);

Dp_5N5 = rho * trapz(yw, uw_5N5 .* (u_5N5 - uw_5N5));    Dp_50 = rho * trapz(yw, uw_50 .* (u_50 - uw_50)); 
Dp_55 = rho * trapz(yw, uw_55 .* (u_55 - uw_55));        Dp_510 = rho * trapz(yw, uw_510 .* (u_510 - uw_510)); 
Dp_515 = rho * trapz(yw, uw_515 .* (u_515 - uw_515));    Dp_520 = rho * trapz(yw, uw_520 .* (u_520 - uw_520)); 

Cdw_5N5 = 2 * trapz(yw,sqrt(qw_5N5/qinf_5N5) - qw_5N5/qinf_5N5);    Cdw_50 = 2 * trapz(yw,sqrt(qw_50/qinf_50) - qw_50/qinf_50); 
Cdw_55 = 2 * trapz(yw,sqrt(qw_55/qinf_55) - qw_55/qinf_55);         Cdw_510 = 2 * trapz(yw,sqrt(qw_510/qinf_510) - qw_510/qinf_510); 
Cdw_515 = 2 * trapz(yw,sqrt(qw_515/qinf_515) - qw_515/qinf_515);    Cdw_520 = 2 * trapz(yw,sqrt(qw_520/qinf_520) - qw_520/qinf_520); 
Cdw = [Cdw_5N5 Cdw_50 Cdw_55 Cdw_510 Cdw_515 Cdw_520];


% xx_5N5 = ; xx_50 = ; 
% xx_55 = ; xx_510 = ; 
% xx_515 = ; xx_520 = ; 

% Theoretical Interpretation
% Initialize 

results = zeros(length(alpha), 5);

for i = 1:length(alpha)
    % Assign temporary variables (_t suffix)
    if alpha(i) == -5
        x_t = N5x; y_t = N5y; Cp_t = N5Cp;
    elseif alpha(i) == 0
        x_t = P0x; y_t = P0y; Cp_t = P0Cp;
    elseif alpha(i) == 5
        x_t = P5x; y_t = P5y; Cp_t = P5Cp;
    elseif alpha(i) == 10
        x_t = P10x; y_t = P10y; Cp_t = P10Cp;
    elseif alpha(i) == 15
        x_t = P15x; y_t = P15y; Cp_t = P15Cp;
    elseif alpha(i) == 20
        x_t = P20x; y_t = P20y; Cp_t = P20Cp;
    end

    % Split into upper and lower surfaces
    upper_mask_t = y_t >= 0;
    x_upper_t = x_t(upper_mask_t);
    y_upper_t = y_t(upper_mask_t);
    Cp_upper_t = Cp_t(upper_mask_t);

    lower_mask_t = ~upper_mask_t;
    x_lower_t = x_t(lower_mask_t);
    y_lower_t = y_t(lower_mask_t);
    Cp_lower_t = Cp_t(lower_mask_t);

    % Ensure ordered from LE to TE (flip upper surface)
    x_upper_t = flip(x_upper_t);
    y_upper_t = flip(y_upper_t);
    Cp_upper_t = flip(Cp_upper_t);

    % Combine into closed loop (upper -> lower)
    x_loop_t = [x_upper_t; x_lower_t];
    y_loop_t = [y_upper_t; y_lower_t];
    Cp_loop_t = [Cp_upper_t; Cp_lower_t];

    % Compute Cn (normal force coefficient)
    dy_t = diff(y_loop_t);
    Cn_t = -sum(Cp_loop_t(1:end-1) .* dy_t);

    % Compute Ca (axial force coefficient)
    dx_t = diff(x_loop_t);
    Ca_t = sum(Cp_loop_t(1:end-1) .* dx_t);

    % Compute Cm_ac (moment about AC)
    x_rel_t = x_loop_t(1:end-1) - xac_c;
    y_rel_t = y_loop_t(1:end-1) - yac_c;
    Cm_ac_t = -sum(Cp_loop_t(1:end-1) .* (x_rel_t .* dy_t - y_rel_t .* dx_t));

    % Compute Cl and Cd (no radians, use cosd/sind)
    alpha_deg_t = alpha(i);
    Cl_t = Cn_t * cosd(alpha_deg_t) - Ca_t * sind(alpha_deg_t);
    Cd_t = Cn_t * sind(alpha_deg_t) + Ca_t * cosd(alpha_deg_t);

    % Store results
    results(i, :) = [alpha_deg_t, Cl_t, Cd_t, Cn_t, Cm_ac_t];
        % results(:,1) results(:,2) results(:,3) results(:,4) results(:,5) 
end

% disp('Results: [Alpha(deg), Cl, Cd, Cn, Cm_ac]');
% disp(results);


% Plots


% AoA vs. Cl
figure(1)
hold on, grid on
plot(alpha, Cl, '-o', "Color", 'r')
plot(results(:,1), results(:,2), 'bo-')
xlabel('Angle of Attack (deg)')
ylabel("C_l")
legend("Experimental", "Theoretical", Location='best')
title("Angle of Attack vs. Lift Coefficient")


% Cl vs. Cd
figure(2)
hold on, grid on
plot(Cl, Cd, '-o', "Color", 'r')
xlabel("C_l")
ylabel("C_d")
title("Lift Coefficient vs. Drag Coefficient")


% Cl vs. Cm_ac
figure(3)
hold on, grid on
plot(Cl, Cm_ac, '-o', "Color", 'r')
xlabel("C_l")
ylabel("Cm_ac")
title("Lift Coefficient vs. Moment Coefficient")


% Cp for AoA = -5
figure(4)
yyaxis left;
hold on, grid on;

hAirfoil = plot(N5x, N5y, 'k-', 'LineWidth', 2); 
fill(N5x, N5y, 'k', 'FaceAlpha', 0.1); 
ylabel('Airfoil Thickness (y/c)');
ylim([-1 1]);

yyaxis right;
hCp = plot(N5x, N5Cp, 'r-', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse'); 
ylabel('Pressure Coefficient (C_p)');
ylim([min(N5Cp)-0.2, max(N5Cp)+0.2]);

plot([Xc_lower' Xc_upper], [Cp_lower_5N5 Cp_upper_5N5], "Color", 'blue', 'LineWidth', 1.5)

xlabel('Chord Position (x/c)');
title('Airfoil Geometry and Pressure Distribution for AoA = -5');
legend('Airfoil', '', "Theoretical", "Experimental", 'Location', 'best');


% Cp for AoA = 0
figure(5)
yyaxis left;
hold on, grid on;

hAirfoil = plot(P0x, P0y, 'k-', 'LineWidth', 2); 
fill(P0x, P0y, 'k', 'FaceAlpha', 0.1); 
ylabel('Airfoil Thickness (y/c)');
ylim([-1 1]);

yyaxis right;
hCp = plot(P0x, P0Cp, 'r-', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse'); 
ylabel('Pressure Coefficient (C_p)');
ylim([min(P0Cp)-0.2, max(P0Cp)+0.2]);

plot([Xc_lower' Xc_upper], [Cp_lower_50 Cp_upper_50], "Color", 'blue', 'LineWidth', 1.5)

xlabel('Chord Position (x/c)');
title('Airfoil Geometry and Pressure Distribution for AoA = 0');
legend('Airfoil', '', "Theoretical", "Experimental", 'Location', 'best');


% Cp for AoA = 5
figure(6)
yyaxis left;
hold on, grid on;

hAirfoil = plot(P5x, P5y, 'k-', 'LineWidth', 2); 
fill(P5x, P5y, 'k', 'FaceAlpha', 0.1); 
ylabel('Airfoil Thickness (y/c)');
ylim([-1 1]);

yyaxis right;
hCp = plot(P5x, P5Cp, 'r-', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse'); 
ylabel('Pressure Coefficient (C_p)');
ylim([min(P5Cp)-0.2, max(P5Cp)+0.2]);

plot([Xc_lower' Xc_upper], [Cp_lower_55 Cp_upper_55], "Color", 'blue', 'LineWidth', 1.5)

xlabel('Chord Position (x/c)');
title('Airfoil Geometry and Pressure Distribution for AoA = 5');
legend('Airfoil', '', "Theoretical", "Experimental", 'Location', 'best');


% Cp for AoA = 10
figure(7)
yyaxis left;
hold on, grid on;

hAirfoil = plot(P10x, P10y, 'k-', 'LineWidth', 2); 
fill(P10x, P10y, 'k', 'FaceAlpha', 0.1); 
ylabel('Airfoil Thickness (y/c)');
ylim([-1 1]);

yyaxis right;
hCp = plot(P10x, P10Cp, 'r-', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse'); 
ylabel('Pressure Coefficient (C_p)');
ylim([min(P10Cp)-0.2, max(P10Cp)+0.2]);

plot([Xc_lower' Xc_upper], [Cp_lower_510 Cp_upper_510], "Color", 'blue', 'LineWidth', 1.5)

xlabel('Chord Position (x/c)');
title('Airfoil Geometry and Pressure Distribution for AoA = 10');
legend('Airfoil', '', "Theoretical", "Experimental", 'Location', 'best');


% Cp for AoA = 15
figure(8)
yyaxis left;
hold on, grid on;

hAirfoil = plot(P15x, P15y, 'k-', 'LineWidth', 2); 
fill(P15x, P15y, 'k', 'FaceAlpha', 0.1); 
ylabel('Airfoil Thickness (y/c)');
ylim([-1 1]);

yyaxis right;
hCp = plot(P15x, P15Cp, 'r-', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse'); 
ylabel('Pressure Coefficient (C_p)');
ylim([min(P15Cp)-0.2, max(P15Cp)+0.2]);

plot([Xc_lower' Xc_upper], [Cp_lower_515 Cp_upper_515], "Color", 'blue', 'LineWidth', 1.5)

xlabel('Chord Position (x/c)');
title('Airfoil Geometry and Pressure Distribution for AoA = 15');
legend('Airfoil', '', "Theoretical", "Experimental", 'Location', 'best');


% Cp for AoA = 20
figure(9)
yyaxis left;
hold on, grid on;

hAirfoil = plot(P20x, P20y, 'k-', 'LineWidth', 2); 
fill(P20x, P20y, 'k', 'FaceAlpha', 0.1); 
ylabel('Airfoil Thickness (y/c)');
ylim([-1 1]);

yyaxis right;
hCp = plot(P20x, P20Cp, 'r-', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse'); 
ylabel('Pressure Coefficient (C_p)');
ylim([min(P20Cp)-0.2, max(P20Cp)+0.2]);

plot([Xc_lower' Xc_upper], [Cp_lower_520 Cp_upper_520], "Color", 'blue', 'LineWidth', 1.5)

xlabel('Chord Position (x/c)');
title('Airfoil Geometry and Pressure Distribution for AoA = 20');
legend('Airfoil', '', "Theoretical", "Experimental", 'Location', 'best');



%Cdw vs. Cdp
figure(10)
hold on, grid on
plot(alpha, Cd, 'Color', 'r')
plot(alpha, Cdw, 'Color', 'b')
xlabel('Angle of Attack (deg)')
ylabel('Drag Coefficient (C_d)')
title('Angle of Attack vs. Drag Coefficient')
legend('Surface Drag', 'Wake Drag', Location='northwest')


% Wake Velocity Profiles
figure(11);
clf; 
hold on;
grid on;

colors = parula(6);  
line_styles = {'-', '--', ':', '-.', '-', '--'};  

legend_labels = {};
q_all = {q5aN5, q5a0, q5a5, q5a10, q5a15, q5a20_1};
alpha_labels = {'-5°', '0°', '5°', '10°', '15°', '20°'};

for i = 1:6
    q_current = q_all{i};
    q_inf = q_current(2) - q_current(1);
    u_wake = sqrt(q_current(33:52)/q_inf); 

    % Plot with style variations
    p = plot(yw, u_wake, ...
        'Color', colors(i,:), ...
        'LineStyle', line_styles{i}, ...
        'MarkerFaceColor', colors(i,:));

    legend_labels{end+1} = ['\alpha = ' alpha_labels{i}];
end

xlabel('Lateral Position, y (m)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Normalized Velocity, u_{wake}/u_{∞}', 'FontSize', 12, 'FontWeight', 'bold');
title('Wake Velocity Profiles vs. Angle of Attack', 'FontSize', 14, 'FontWeight', 'bold');

ylim([0.4 1.1]);
xlim([min(yw) max(yw)]);
legend(legend_labels)