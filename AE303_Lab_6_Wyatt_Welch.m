%% Description -- Please Read!
%
% This script imports the data from the .csv files for the PIV lab of
% AE303. A status is printed to the command window as each file is imported
% and then once again as each file's data is formatted. Each .csv file
% represents one instant in time for the flow. The data was recorded at 150
% Hz, so 300 files corresponds to 2 seconds of data. 
%
% This script checks if the data has been read before so that one does not
% spend extraneous time re-importing data. However, if the import portion
% of the script is terminated early, the parsing portion of the script will
% not function properly and the workspace should be cleared before running
% this script again. 
%
% To further avoid repeating portions of the script, divide the script
% into sections with the '%% Title' as is done already and take advantage
% of the 'Run Section' option to the right of 'Run' in Matlab.
%
% The purpose of this script is to allow the student to work with the data
% presented rather than put unnecessary time into writing their own import
% and formatting script since there is such a large amount of data.
%
%% Notes on the script output
%
% Upon this script's completion, there will exist 10 new variables in the
% workspace. Each is important for processing the PIV data. They are as
% follows:
%
%   data        This is the raw 3D array containing all of the data from
%               all 300 of the .csv files; dimensions 8806x11x300
%
%   N           The number of time steps in the dataset; 300
%
%   winsize     The 1x2 array with the number of rows and columns,
%               respectively, in the vector grid; should be [119,74]
%
%   X           The 2D matrix of the x-coordinates of the vector grid
%
%   Y           The 2D matris of the y-coordinates of the vector grid
%
%   u           The 2D matrix of the u-component of the velocity on the
%               vector grid with likely erroneous vectors replaced by NaN
%               with units of m/s
%
%   v           The 2D matrix of the v-component of the velocity on the
%               vector grid with likely erroneous vectors replaced by NaN
%               with units of m/s
%
%   u_raw       The 2D matrix of the u-component of the velocity on the
%               vector grid with no replaced vectors and units of m/s
%
%   v_raw       The 2D matrix of the v-component of the velocity on the
%               vector grid with no replaced vectors and units of m/s
%
%   u_pix       The 2D matrix of the u-component of the velocity on the
%               vector grid with likely erroneous vectors replaced by NaN
%               with units of pixel displacement
%
%   v_pix       The 2D matrix of the v-component of the velocity on the
%               vector grid with likely erroneous vectors replaced by NaN
%               with units of pixel displacement
%
% The velocity output arrays are all 3D arrays, where the first two
% dimensions represent the vector grid and the third dimension is how each
% time step is stored. Moving, for example, from u(:,:,1) to u(:,:,2) is 
% moving one time step forward in time, while referencing the full vector
% grid of u velocity components.
%
% The bottom of the parse portion of the script includes a commented out
% code which will plot and animation of the flow as the data is parsed.
% Running this portion can cause the code to take longer than desired to
% run, which is why it is commented out.
%
%% Read Data
if ~exist('data','var') % Check if data read already

    winsize = [119 74]; % Size of vector grid
    L = winsize(1)*winsize(2); % Get size of vector field

    N = 300; % Number of files (time steps)
    
    data = zeros(L,11,N); % Pre-allocate

    for i = 1:N
        fprintf('Now Retrieving File %d of 300\n',i);
        file = sprintf('EduPIV_lab.62tbxosb.%06d.csv',i-1);
        data(:,:,i) = readmatrix(file); % Save each file to data
    end
end

%% Parse Data
s = 1;
t = 1;
u = zeros(119,74,300); v = u;
u_raw = u; u_pix = u;
v_raw = v; v_pix = v; % Pre-Allocate
x = u; y = v;
N = 300;
for i = 1:N
    fprintf('Now Parsing Set %d of 300\n',i); % Status display
    for q = 1:8806
        if data(q,11,i) ~= 0
            u(s,t,i) = nan; % Replace flagged vectors with nan
            v(s,t,i) = nan;
            u_pix(s,t,i) = nan;
            v_pix(s,t,i) = nan;
        else
            u(s,t,i) = data(q,9,i); % Parse valid vectors
            v(s,t,i) = data(q,10,i);
            u_pix(s,t,i) = data(q,7,i);
            v_pix(s,t,i) = data(q,8,i);
        end
        
        u_raw(s,t,i) = data(q,9,i);
        v_raw(s,t,i) = data(q,10,i);
        
        x(s,t,i) = data(q,5,i); % Record grid
        y(s,t,i) = data(q,6,i);
        
        s = s + 1; % Move to next row
        
        if mod(q,winsize(1)) == 0 % Detect is done with column
            t = t + 1; % Move to next vector column
            s = 1; % Reset to first row
        end
        if q == L
            t = 1; % Finished with frame, reset t
%             figure(1) % Update only figure 1
%             contourf(x(:,:,i),y(:,:,i),sqrt(u(:,:,i).^2 + v(:,:,i).^2),20,'LineStyle','none')
%             axis equal % Update cool animation of the flow
%             xlabel('x [cm]')
%             ylabel('y [cm]')
%             drawnow
        end
    end
end

X = x(:,:,1); % Grid doesn't change, so we just need a 2D matrix
Y = y(:,:,1);

clear file i q s t x y L

%% Data Processing
% Calculation
D = 0.05; % m
mu = 0.0010016;
rho = 998.23;
Xndm = X / D;
Yndm = Y / D;

uAvg = mean(u, 3, 'omitnan');
vAvg = mean(v, 3, 'omitnan');

uDif = u - uAvg;
vDif = v - vAvg;

uStr = mean(uDif .^ 2, 3, 'omitnan');
vStr = mean(vDif .^ 2, 3, 'omitnan');

uvStr = mean(uDif .* vDif, 3, 'omitnan');

[dudy, dudx] = gradient(uAvg, Y(1,:), X(:,1));
[dvdy, dvdx] = gradient(vAvg, Y(1,:), X(:,1));

vort = dvdx - dudy;

uex = uAvg(1, 12:59);
vex = vAvg(1, 12:59);

Ucenter = mean([uAvg(:,35), uAvg(:,36)], 'all');
Vcenter = mean([vAvg(:,35), vAvg(:,36)], 'all');

Uex = mean(sqrt(uex .^ 2 + vex .^ 2));

Re = rho * Uex * D / mu;

Ycenter = Y / D - (max(Y / D) - min(Y / D)) / 2; % Centered

%% Plots

figure(1)
hold on, grid on
contourf(Xndm, Yndm, uAvg, 50, 'LineColor', 'none')

colorbar;
ylabel(colorbar, "Time Averaged u-Velocity (m/s)")
xlabel('X / Diameter')
ylabel('Y / Diameter')
title('Time Averaged u-Velocity Contour Map')



figure(2)
hold on, grid on
contourf(Xndm, Yndm, uStr, 50, 'LineColor', 'none')

colorbar;
ylabel(colorbar, "Reynolds Stress u'² (m²/s²)")
xlabel('X / Diameter')
ylabel('Y / Diameter')
title('Reynolds Stress $\overline{u''^2}$', 'Interpreter', 'latex')



figure(3)
hold on, grid on
contourf(Xndm, Yndm, vStr, 50, 'LineColor', 'none')

colorbar;
ylabel(colorbar, "Reynolds Stress v'² (m²/s²)")
xlabel('X / Diameter')
ylabel('Y / Diameter')
title('Reynolds Stress $\overline{v''^2}$', 'Interpreter', 'latex')



figure(4)
hold on, grid on
contourf(Xndm, Yndm, uvStr, 50, 'LineColor', 'none')
ylabel(colorbar, "Reynolds Stress u'²v'² (m²/s²)")

colorbar;
ylabel(colorbar, "Reynolds Shear Stress u'²v'² (m²/s²)")
xlabel('X / Diameter')
ylabel('Y / Diameter')
title('Reynolds Shear Stress $\overline{u''v''^2}$', 'Interpreter', 'latex')



figure(5)
hold on, grid on
contourf(Xndm, Yndm, vort, 50, 'LineColor', 'none')

colorbar;
ylabel(colorbar, "Vorticity Magnitude $1/s$", 'Interpreter', 'latex')
xlabel('X / Diameter')
ylabel('Y / Diameter')
title('Vorticity Magnitude Disruption')



figure(6)
hold on, grid on
plot(Ycenter(1,:), uAvg(10,:) / Ucenter, 'r')
plot(Ycenter(1,:), uAvg(40,:) / Ucenter, 'g')
plot(Ycenter(1,:), uAvg(80,:) / Ucenter, 'b')
plot(Ycenter(1,:), uAvg(119,:) / Ucenter, 'm')

xlabel("Y / Diameter")
ylabel("Velocity / Centerline of u-Velocity")
title("Streamwise Variations of u-Component Velocity")
legend("$X/D = 0.2$", "$X/D = 0.9$", "$X/D = 1.8$", "$X/D = 2.7$", "Interpreter", "latex");



figure(7)
hold on, grid on
plot(Ycenter(1,:), vAvg(10,:) / Ucenter, 'r')
plot(Ycenter(1,:), vAvg(40,:) / Ucenter, 'g')
plot(Ycenter(1,:), vAvg(80,:) / Ucenter, 'b')
plot(Ycenter(1,:), vAvg(119,:) / Ucenter, 'm')

xlabel("Y / Diameter")
ylabel("Velocity / Centerline of v-Velocity")
title("Streamwise Variations of v-Component Velocity")
legend("$X/D = 0.2$", "$X/D = 0.9$", "$X/D = 1.8$", "$X/D = 2.7$", "Interpreter", "latex");



figure(8)
hold on, grid on
plot(Ycenter(1,:), uStr(10,:), 'r')
plot(Ycenter(1,:), uStr(40,:), 'g')
plot(Ycenter(1,:), uStr(80,:), 'b')
plot(Ycenter(1,:), uStr(119,:), 'm')

xlabel("Y / Diameter")
ylabel("Reynolds Stress u'² (m²/s²)")
title("Streamwise Variations of u-Component Reynolds Stress")
legend("$X/D = 0.2$", "$X/D = 0.9$", "$X/D = 1.8$", "$X/D = 2.7$", "Interpreter", "latex");



figure(9)
hold on, grid on
plot(Ycenter(1,:), vStr(10,:), 'r')
plot(Ycenter(1,:), vStr(40,:), 'g')
plot(Ycenter(1,:), vStr(80,:), 'b')
plot(Ycenter(1,:), vStr(119,:), 'm')

xlabel("Y / Diameter")
ylabel("Reynolds Stress v'² (m²/s²)")
title("Streamwise Variations of v-Component Reynolds Stress")
legend("$X/D = 0.2$", "$X/D = 0.9$", "$X/D = 1.8$", "$X/D = 2.7$", "Interpreter", "latex");



figure(10)
hold on, grid on
plot(Ycenter(1,:), uvStr(10,:), 'r')
plot(Ycenter(1,:), uvStr(40,:), 'g')
plot(Ycenter(1,:), uvStr(80,:), 'b')
plot(Ycenter(1,:), uvStr(119,:), 'm')

xlabel("Y / Diameter")
ylabel("Reynolds Shear Stress u'²v'² (m²/s²)")
title("Streamwise Variations of Reynolds Shear Stress")
legend("$X/D = 0.2$", "$X/D = 0.9$", "$X/D = 1.8$", "$X/D = 2.7$", "Interpreter", "latex", 'Location', 'northwest');



figure(11)
hold on, grid on
plot(Ycenter(1,:), vort(10,:), 'r')
plot(Ycenter(1,:), vort(40,:), 'g')
plot(Ycenter(1,:), vort(80,:), 'b')
plot(Ycenter(1,:), vort(119,:), 'm')

xlabel("Y / Diameter")
ylabel("Vorticity Magnitude $1/s$", 'Interpreter', 'latex')
title("Streamwise Variations of Vorticity")
legend("$X/D = 0.2$", "$X/D = 0.9$", "$X/D = 1.8$", "$X/D = 2.7$", "Interpreter", "latex", 'Location', 'northwest');