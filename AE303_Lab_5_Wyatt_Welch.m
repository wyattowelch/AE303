%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AE302 Lab 5 - Wyatt Welch
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear all

% Load Data
load("dat111.mat"); dat111 = dat111(2:end,:);
load("dat110.mat"); 
load("dat101.mat"); dat101 = dat101(2:end,:);
load("dat001.mat"); 
load("dat000.mat"); 

% Corrections 
fill001 = dat001(1,:);
dat001 = [repmat(fill001, 10, 1); dat001(2:3,:)];
fill000 = dat000(1,:);
dat000 = [repmat(fill000, 10, 1); dat000(2:3,:)];

% Reordered Data
q = 5 * .0360912; %Psi
S = 93.81;
c = 3.466;
b = 27.066;
R = 1716; %ftlb/slugR
angA = dat111(:,1);
angB = dat000(10:12,2);

T111 = 75.1; P111 = 30.06;
T110 = 79; P110 = 29.97;
T101 = 72.5; P101 = 30.08;
T001 = 78.8; P001 = 29.96;
T000 = 79.1; P000 = 29.97;
T_array = [T111 T110 T101 T001 T000] + 459.67; %R
P_array = [P111 P110 P101 P001 P000] .* .49115; %Psi

% Reynolds Numbers
rho = (P_array * 144) ./ (R * T_array);
V = sqrt((2 * (5 * 5.204)) ./ rho);
mu = 2.629e-7 * (T_array.^1.5) ./ (T_array + 198.72);
Re = (rho .* V * (c / 12)) ./ mu;
fprintf("Reynolds numbers [111 110 101 001 000]:\n")
fprintf("%4.4e ",Re)
fprintf("\n\n")

% Set-up Data
Fxd111 = dat111(:,4);
Fxd110 = dat110(:,4);
Fxd101 = dat101(:,4);
Fxd001 = dat001(:,4);
Fxd000 = dat000(:,4);

Fyd111 = dat111(:,5);
Fyd110 = dat110(:,5);
Fyd101 = dat101(:,5);
Fyd001 = dat001(:,5);
Fyd000 = dat000(:,5);

Fzd111 = dat111(:,6);
Fzd110 = dat110(:,6);
Fzd101 = dat101(:,6);
Fzd001 = dat001(:,6);
Fzd000 = dat000(:,6);

Mxd111 = dat111(:,7);
Mxd110 = dat110(:,7);
Mxd101 = dat101(:,7);
Mxd001 = dat001(:,7);
Mxd000 = dat000(:,7);

Myd111 = dat111(:,8);
Myd110 = dat110(:,8);
Myd101 = dat101(:,8);
Myd001 = dat001(:,8);
Myd000 = dat000(:,8);

Mzd111 = dat111(:,9);
Mzd110 = dat110(:,9);
Mzd101 = dat101(:,9);
Mzd001 = dat001(:,9);
Mzd000 = dat000(:,9);

% Force/Moment Corrections
FxT1 = (Fxd111 - Fxd110) - (Fxd001 - Fxd000);
FyT1 = (Fyd111 - Fyd110) - (Fyd001 - Fyd000);
FzT1 = (Fzd111 - Fzd110) - (Fzd001 - Fzd000);

MxT1 = (Mxd111 - Mxd110) - (Mxd001 - Mxd000);
MyT1 = (Myd111 - Myd110) - (Myd001 - Myd000);
MzT1 = (Mzd111 - Mzd110) - (Mzd001 - Mzd000);


FxT0 = (Fxd101 - Fxd110) - (Fxd001 - Fxd000);
FyT0 = (Fyd101 - Fyd110) - (Fyd001 - Fyd000);
FzT0 = (Fzd101 - Fzd110) - (Fzd001 - Fzd000);

MxT0 = (Mxd101 - Mxd110) - (Mxd001 - Mxd000);
MyT0 = (Myd101 - Myd110) - (Myd001 - Myd000);
MzT0 = (Mzd101 - Mzd110) - (Mzd001 - Mzd000);


% Coefficient Calculations
Cl1 = FzT1 / (q * S);
Cl0 = FzT0 / (q * S);

Cd1 = FxT1 / (q * S);
Cd0 = FxT0 / (q * S);

Cm1 = MyT1 / (q * S * c);
Cm0 = MyT0 / (q * S * c);

Cn1 = MzT1 / (q * S * b);
Cn0 = MzT0 / (q * S * b);

% Cd Components
AR = b^2 / S;
e = 1.78 * (1-.045 * AR ^ .68) - .64;
fprintf("The Oswald Efficiency Factor, e = %4.4f\n\n",e)
K = 1/(pi() * e * AR);

alpha_dat = angA(1:end-2);
Cl_dat = Cl1(1:end-2);
Cl_dat0 = Cl0(1:end-2);
cross = find(diff(sign(Cl_dat)));
cross0 = find(diff(sign(Cl_dat0)));
alpha0 = interp1(Cl_dat(cross:cross+1), ...
               alpha_dat(cross:cross+1), ...
               0, 'linear');
alpha01 = interp1(Cl_dat0(cross0:cross0+1), ...
               alpha_dat(cross0:cross0+1), ...
               0, 'linear');
fprintf("Zero lift angles \nTail On: %4.4f\nTail Off: %4.4f\n\n",alpha0,alpha01)
[~,idx] = min(abs(angA-alpha0));
C_d0 = Cd1(idx);
Cd = C_d0 + K * Cd1.^2;

ClCdMax1 = max(max(Cl1(1:end-2) / Cd1(1:end-2)));
ClCdMax0 = max(max(Cl0(1:end-2) / Cd0(1:end-2)));
fprintf("Cl/Cd maximums\nTail on: %4.4f \nTail off: %4.4f\n\n",ClCdMax1,ClCdMax0)

% dCl/da, dCm/da, dCn/db
lin_range = (alpha_dat >= -6) & (alpha_dat <=8);
pCl = polyfit(alpha_dat(lin_range), Cl_dat(lin_range),1);
Clslope = pCl(1);

Cm_dat = Cm1(1:end-2);
pCm = polyfit(alpha_dat(lin_range), Cm_dat(lin_range),1);
Cmslope = pCm(1);

beta_dat = angB;
lin_range = (beta_dat >= 0) & (beta_dat <=10);
Cn_dat = Cn1(10:12,:);
pCn = polyfit(beta_dat(lin_range), Cn_dat(lin_range),1);
Cnslope = pCn(1);

fprintf("Lift Slope dCl/da = %4.4f \nPitch Moment Slope dCm/da = %4.4f " + ...
    "\nYaw Moment Slope dCn/db = %4.4f\n\n",Clslope,Cmslope, Cnslope)


% Stall
[Cl_max, idx_max] = max(Cl_dat);
alpha_stall = alpha_dat(idx_max);
fprintf("Maximum lift coefficient, C_lmax = %4.4f" + ...
    "\nCritical angle of attack, \alpha_max = %4.4f\n\n",Cl_max,alpha_stall)


% Plots
figure(1)
hold on, grid on
plot(angA(1:end-2), Cl1(1:end-2), '-ob')
plot(angA(1:end-2), Cl0(1:end-2), '-or')

legend("Tail On", "Tail Off", 'Location','northwest')
xlabel('Angle of Attack (\alpha)')
ylabel('Lift Coefficient (C_L)')
title('Angle of Attack vs. Lift Coefficient')


figure(2)
hold on, grid on
plot(angA(1:end-2), Cm1(1:end-2), '-ob')
plot(angA(1:end-2), Cm0(1:end-2), '-or')

legend("Tail On", "Tail Off", 'Location','northeast')
xlabel('Angle of Attack (\alpha)')
ylabel('Pitching Moment Coefficient (C_M)')
title('Angle of Attack vs. Pitching Moment Coefficient')


figure(3)
hold on, grid on
plot(angB, Cn1(10:12,:), '-ob')
plot(angB, Cn0(10:12,:), '-or')

legend("Tail On", "Tail Off", 'Location','northeast')
xlabel('Sideslip Angle (\beta)')
ylabel('Yawing Moment Coefficient (C_N)')
title('Sideslip Angle vs. Yawing Moment Coefficient')


figure(4)
hold on, grid on
plot(Cd1(1:end-2), Cl1(1:end-2), '-ob')
plot(Cd0(1:end-2), Cl0(1:end-2), '-or')

legend("Tail On", "Tail Off", 'Location','northwest')
ylabel('Lift Coefficient (C_L)')
xlabel('Drag Coefficient (C_D)')
title('Lift Coefficient vs. Drag Coefficient')