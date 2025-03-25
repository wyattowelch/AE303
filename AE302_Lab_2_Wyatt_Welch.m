%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AE302 Lab 2 - Wyatt Welch
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear all, load('RawArrays.mat')

S1 = table2array(S1); % PSI
S2 = table2array(S2); 
S3 = table2array(S3);

% Remove last 2 columns, zero col 7, replace col 2 with average
A1 = S1(2:end, 1:end-2);
A2 = S2(2:end, 1:end-2);
A3 = S3(2:end, 1:end-2);

B1 = A1 - A1(:,7);
B2 = A2 - A2(:,7);
B3 = A3 - A3(:,7);

C1 = mean(B1, 2);
C2 = mean(B2, 2);
C3 = mean(B3, 2);

B1(:,2) = C1;
B2(:,2) = C2;
B3(:,2) = C3;

S1c = mean(B1, 1);
S2c = mean(B2, 1);
S3c = mean(B3, 1);

% Data
AmT0 = 71.7 + 459.67; % F to R
WtT0 = 78.4 + 459.67;
AmP0 = 30.16 / 2.036; % inHg to PSI

AmT2 = 71.7 + 459.67;
WtT2 = 76.6 + 459.67;
AmP2 = AmP0;

AmT5 = 71.7 + 459.67;
WtT5 = 82.7 + 459.67;
AmP5 = AmP0;

% Time-Averaging Dynamic Pressure for Each Pitot Tube

S1Avg=mean(B1,1);
S2Avg=mean(B2,1);
S3Avg=mean(B3,1);

q_hat1 = mean(S1Avg);
q_hat2 = mean(S2Avg);
q_hat3 = mean(S3Avg);

% Average pressure at given dynamic pressure setting
S1_total = B1;
S2_total = B2;
S3_total = B3;

P_avg1 = mean(S1_total, 2);  % Average across columns (all pitots)
P_avg2 = mean(S2_total, 2);
P_avg3 = mean(S3_total, 2);

Static1 = S1_total(:,7);
Static2 = S2_total(:,7);
Static3 = S3_total(:,7);

q_avg1 = mean(P_avg1 - Static1);
q_avg2 = mean(P_avg2 - Static2);
q_avg3 = mean(P_avg3 - Static3);

% Calculate free stream velocity, Reynold's Number

den1 = (AmP0) / (1716 * WtT0); % slugs/ft³
den2 = (AmP2) / (1716 * WtT2);
den3 = (AmP5) / (1716 * WtT5);

V_avg1 = sqrt((2*S1Avg)/den1);
V_avg2 = sqrt((2*S2Avg)/den2);
V_avg3 = sqrt((2*S3Avg)/den3);

V_hat1 = sqrt((2*mean(S1c([1:6, 8:end]), 2))/den1);
V_hat2 = sqrt((2*mean(S2c([1:6, 8:end]), 2))/den2);
V_hat3 = sqrt((2*mean(S3c([1:6, 8:end]), 2))/den3);

V_theo1 = sqrt((2*0)/den1);
V_theo2 = sqrt((2*0.0721824/den2));
V_theo3 = sqrt((2*0.180456/den3));

mu_ref = 3.737e-7; % Reference dynamic viscosity in slugs/(ft·s)
T_ref = 518.67;    % Reference temperature in Rankine
S = 198.72;        % Sutherland's constant in Rankine

mu0 = mu_ref * (WtT0 / T_ref)^(3/2) * ((T_ref + S) / (WtT0 + S));
mu2 = mu_ref * (WtT2 / T_ref)^(3/2) * ((T_ref + S) / (WtT2 + S));
mu5 = mu_ref * (WtT5 / T_ref)^(3/2) * ((T_ref + S) / (WtT5 + S));

Re1 = (den1 * V_hat1) / mu0;
Re2 = (den2 * V_hat2) / mu2;
Re3 = (den3 * V_hat3) / mu5;

% Plots
pitotN = 1:35;
N = size(pitotN);


figure(1) %Flow uniformity - Measured Pressure

hold on, grid on
ylim([-.02 .2])

scatter(pitotN,S1c)
plot(pitotN, ones(1, length(pitotN)) * mean(S1c(:, [1:6, 8:end]), 2), 'k')
scatter(pitotN,S2c)
plot(pitotN, ones(1, length(pitotN)) * mean(S2c(:, [1:6, 8:end]), 2), 'k')
scatter(pitotN,S3c)
plot(pitotN, ones(1, length(pitotN)) * mean(S3c(:, [1:6, 8:end]), 2), 'k')

xlabel("Pressure Ports")
ylabel("Pressure (Psi)")
title("Flow Uniformity - Measured Pressure")
legend("0.0 inH2O", '', "2.0 inH2O", '', "5.0 inH2O", "Average")


figure(2) %Flow uniformity - Dynamic Pressure
hold on, grid on
ylim([-.02 .22])

scatter(pitotN, S1Avg)
plot(pitotN, ones(N) * mean(S1c(:, [1:6, 8:end]), 2), 'r')
plot(pitotN, ones(N) * 0, 'k')
scatter(pitotN, S2Avg)
plot(pitotN, ones(N) * mean(S2c(:, [1:6, 8:end]), 2), 'r')
plot(pitotN, ones(N) * 0.0721825, 'k')
scatter(pitotN, S3Avg)
plot(pitotN, ones(N) * mean(S3c(:, [1:6, 8:end]), 2), 'r')
plot(pitotN, ones(N) * .180456, 'k')

xlabel("Pressure Ports")
ylabel("Pressure (Psi)")
title("Flow Uniformity - Dynamic Pressure")
legend("0.0 inH2O", '', '', "2.0 inH2O", '', '', "5.0 inH2O", "Average", 'q_s_e_t_t_i_n_g')


figure(3)
hold on, grid on
ylim([-20 170])

scatter(pitotN,V_avg1)
plot(pitotN, ones(N) * V_hat1,'r')
plot(pitotN, ones(N) * V_theo1, 'k')
scatter(pitotN,V_avg2)
plot(pitotN, ones(N) * V_hat2,'r')
plot(pitotN, ones(N) * V_theo2, 'k')
scatter(pitotN,V_avg3)
plot(pitotN, ones(N) * V_hat3,'r')
plot(pitotN, ones(N) * V_theo3, 'k')


xlabel("Pressure Ports")
ylabel("Airspeed (ft/s)")
title("Flow Uniformity - Test Section Airspeed")
legend("0.0 inH2O", '', '', "2.0 inH2O", '', '', "5.0 inH2O", "Average", "v_s_e_t_t_i_n_g")


figure(4)
hold on, grid on
ylim([-1 6])

scatter(pitotN, 100 * (S2Avg - mean(S2Avg)) / mean(S2Avg),'k');
scatter(pitotN, 100 * (S3Avg - mean(S3Avg)) / mean(S3Avg),'r');

xlabel("Pressure Ports")
ylabel('$\frac{(q - \bar{q})}{\bar{q}}$', 'Interpreter', 'latex');
title("Flow Uniformity - Dynamic Pressure Deviation From Average")
legend("2.0 inH2O", "5.0 inH2O")
ytickformat('percentage')


figure(5)
hold on, grid on
ylim([-1 6])

scatter(pitotN, 100 * (S2Avg - mean(S2Avg)) / mean(S2Avg),'k');

xlabel("Pressure Ports")
ylabel('$\frac{(q - \bar{q})}{\bar{q}}$', 'Interpreter', 'latex');
title("Flow Uniformity - Dynamic Pressure Deviation From Average")
legend("2.0 inH2O")
ytickformat('percentage')


figure(6)
hold on, grid on
ylim([-1 6])

scatter(pitotN, 100 * (S3Avg - mean(S3Avg)) / mean(S3Avg),'r');

xlabel("Pressure Ports")
ylabel('$\frac{(q - \bar{q})}{\bar{q}}$', 'Interpreter', 'latex');
title("Flow Uniformity - Dynamic Pressure Deviation From Average")
legend("5.0 inH2O")
ytickformat('percentage')


figure(7)
hold on, grid on
ylim([1 5])

scatter(pitotN, 100 * (V_avg2 - mean(V_avg2)) / mean(V_avg2), 'k')
scatter(pitotN, 100 * (V_avg3 - mean(V_avg3)) / mean(V_avg3), 'r')

xlabel("Pressure Ports")
ylabel('$\frac{(V - \bar{V})}{\bar{V}}$', 'Interpreter', 'latex');
title("Flow Uniformity - Airspeed Deviation From Average")
legend("2.0 inH2O", "5.0 inH2O")
ytickformat('percentage')


figure(8)
hold on, grid on
ylim([1 5])

scatter(pitotN, 100 * (V_avg2 - mean(V_avg2)) / mean(V_avg2), 'k')

xlabel("Pressure Ports")
ylabel('$\frac{(V - \bar{V})}{\bar{V}}$', 'Interpreter', 'latex');
title("Flow Uniformity - Airspeed Deviation From Average")
legend("2.0 inH2O")
ytickformat('percentage')


figure(9)
hold on, grid on
ylim([1 5])

scatter(pitotN, 100 * (V_avg3 - mean(V_avg3)) / mean(V_avg3), 'r')

xlabel("Pressure Ports")
ylabel('$\frac{(V - \bar{V})}{\bar{V}}$', 'Interpreter', 'latex');
title("Flow Uniformity - Airspeed Deviation From Average")
legend("5.0 inH2O")
ytickformat('percentage')

%Dynamic Pressure Deviation from Average (q = 2.0 inH2O)
x = [linspace(-19.5, 19.5, 15), repelem(-7.5,10), repelem(7.5,10)]';
z = [repelem(0,15), (-13.5:3:13.5), (-13.5:3:13.5)]';
S2Avg = S2Avg';
[X,Z] = meshgrid(linspace(min(x), max(x), 35), linspace(min(z), max(z),35));
F2 = scatteredInterpolant(x,z, (S2Avg - mean(S2Avg)) / mean(S2Avg), "natural", "none");
QD2 = F2(X,Z);
QD2(QD2 < 0) = .01; 
figure(10)
hold on, grid on
contourf(X,Z,QD2,20,'LineColor','none')
title('Dynamic Pressure Deviation from Average (q = 2.0 inH2O) (%)')
xlabel("X Position (in)")
ylabel("Z Position (in)")
cb = colorbar;
ylabel(cb, '$\frac{(q - \bar{q})}{\bar{q}}$', 'Interpreter', 'latex');
cb.TickLabels = strcat(cb.TickLabels, '%');
plot(x,z,'ro')
xlim([-20 20])
ylim([-15 15])
clim([0 .06])

S3Avg = S3Avg';
F3 = scatteredInterpolant(x,z, (S3Avg - mean(S3Avg)) / mean(S3Avg), "natural", "none");
QD3 = F3(X,Z);
QD3(QD3 < 0) = .01; 
figure(11)
hold on, grid on
contourf(X,Z,QD3,20,'LineColor','none')
title('Dynamic Pressure Deviation from Average (q = 5.0 inH2O) (%)')
xlabel("X Position (in)")
ylabel("Z Position (in)")
cb = colorbar;
ylabel(cb, '$\frac{(q - \bar{q})}{\bar{q}}$', 'Interpreter', 'latex');
cb.TickLabels = strcat(cb.TickLabels, '%');
plot(x,z,'ro')
xlim([-20 20])
ylim([-15 15])
clim([-.02 .06])


IAS = 54.5 * 3.28084;
V_IAS = sqrt((2*S3Avg)/den3);
V_IAS(V_IAS <1) = 153.2294;
F4 = scatteredInterpolant(x,z, V_IAS, "natural", "none");
QD4 = F4(X,Z);
QD3(QD3 < 1) = 153.2294; 
figure(12)
hold on, grid on
contourf(X,Z,QD4,20,'LineColor','none')
title('Airspeed Distribution - IAS = 54.5 m/s = 178.8 ft/s')
xlabel("X Position (in)")
ylabel("Z Position (in)")
cb = colorbar;
ylabel(cb, 'Measured Airspeed (ft/s)');
plot(x,z,'ro')
xlim([-20 20])
ylim([-15 15])
clim([150 156])


% Convergence
B17 = B1;
B17(:,7) = [];
N = length(B17);
threshold = 0.00005;
cum_mean = zeros(N,1);
perc_conv = zeros(N,1);

for i = 1:N
    cum_mean(i) = mean(B17(1:i));
    if i == 1
        perc_conv(1) = NaN;
    else
        perc_conv(i) = abs((cum_mean(i) - cum_mean(i-1)) / cum_mean(i-1));
    end
end

conv_idx = find(perc_conv < threshold, 1);
if isempty(conv_idx)
    conv_idx = N;
    fprintf('not converged')
else
    fprintf('converged at simple %d with %% change = %f\n', conv_idx,perc_conv(conv_idx));
end

figure(13)
subplot(2,1,1);
hold on, grid on
plot(1:N, cum_mean)
plot(conv_idx, cum_mean(conv_idx))
xlabel('Number of Samples')
ylabel('Cumulative Mean Dynamic Pressure')
title('Convergence Test for 0.0 inH2O')
legend('Cumulative Mean', sprintf('Convergence at n=%d',conv_idx))

subplot(2,1,2)
plot(2:N, perc_conv(2:end))
hold on
plot(conv_idx, perc_conv(conv_idx))
xlabel('Number of Samples')
ylabel('Percentage Convergence')
title('Percentage Convergence Criterion')
yline(threshold, 'r--')
legend('Percentage Convergence', sprintf('Threshold = %.5f', threshold), "Threshold")


B27 = B2;
B27(:,7) = [];
N = length(B27);
threshold = 0.00005;
cum_mean = zeros(N,1);
perc_conv = zeros(N,1);

for i = 1:N
    cum_mean(i) = mean(B27(1:i));
    if i == 1
        perc_conv(1) = NaN;
    else
        perc_conv(i) = abs((cum_mean(i) - cum_mean(i-1)) / cum_mean(i-1));
    end
end

conv_idx = find(perc_conv < threshold, 1);
if isempty(conv_idx)
    conv_idx = N;
    fprintf('not converged')
else
    fprintf('converged at simple %d with %% change = %f\n', conv_idx, perc_conv(conv_idx));
end

figure(14)
subplot(2,1,1);
hold on, grid on
plot(1:N, cum_mean)
plot(conv_idx, cum_mean(conv_idx))
xlabel('Number of Samples')
ylabel('Cumulative Mean Dynamic Pressure')
title('Convergence Test for 2.0 inH2O')
legend('Cumulative Mean', sprintf('Convergence at n=%d', conv_idx))

subplot(2,1,2)
plot(2:N, perc_conv(2:end))
hold on
plot(conv_idx, perc_conv(conv_idx))
xlabel('Number of Samples')
ylabel('Percentage Convergence')
title('Percentage Convergence Criterion')
yline(threshold, 'r--')
legend('Percentage Convergence', sprintf('Threshold = %.5f', threshold), "Threshold")
