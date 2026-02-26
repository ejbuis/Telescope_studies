% Script to check the decay of a spherical and cylindrical source against
% expected literature values as well as the decay of a typical shower

clear all; 
addpath('..');
addpath('../Auxiliary');
addpath('Showerparam_extensions');
addpath('Datasets');
rng(10); % FIX RANDOM SEED!

Eo = 1e9;           % Primary Energy        [GeV],  only for Corsika
plot_bool = 1;      % Whether to plot       [-]
save_bool = 0;      % Whether to save       [-]
nmc = 1e6;          % Number of MC points   [-],    default 1e6
fs = 1e6;           % Sampling frequency    [Hz],   default 1e6 Hz
c = 1500;           % speed of sound        [m/s]
wind = 4096*4;      % 1/2 size of the time window [-]
atten = 5;          % Learned's attenuation = 1, No attenuation = 5
debug = 1;

% Source specifics
z_max = 0; 	                % Energy Distribution z coordinate OR length of cylinder
rpos = logspace(1,5,20);    % Observer r position   [m],    default 1000 m
zpos = 6;                   % Observer z position   [m],    default 6 m

sigma = 1;                  % Energy Distribution Maximum Standard deviation [cm]
w = 0;                      % Weight maximum vs cilinder [a.u.]

tic

%% ============== Generate profile from shower ==============
t_axis = (-wind:wind-1)/fs; % time axis (centered at 0)

% => Binning specifics
rsc = [.5:9.5 15:10:105];   % Radial bin centres (cm)
zsc = 10:20:2000;            % Longitudinal Bin Centres (cm)

tsmc = ShowerParm(rsc, zsc, Eo, 'CORSIKA');

% as the 10-100cm bins are 10x wider need to scale by a factor of 10 
tsmc = tsmc * diag(kron([1 10], ones(1, 10)));
pointsc = MCGEn(tsmc, [0 zsc+10], [0 rsc+[0.5*ones(1, 10) 5*ones(1, 10)]],nmc);

% Convert to cartesian
[x,y,z] = pol2cart(rand(nmc, 1)*2*pi, pointsc(:, 2), pointsc(:, 1));

points = [x y z] * 1e-2; 
pmaxs = zeros(length(rpos),1);

for i=1:length(rpos)
    Do  = [rpos(i) 0 zpos];        % Position of observer
    [tplot, p, pw, E] = kernelfrmod(points, Do, log10(Eo), atten, 5, fs, wind);
    pmaxs(i) = max(p);
end

%% Fit appropriate curve
f1 = fit(rpos', pmaxs,'a/(x-b) + c/sqrt(x-d)', 'StartPoint', [1,1,1,1]);
r_fit1 = logspace(log10(min(rpos)), log10(max(rpos)), 100);
y_fit1 = f1(r_fit1);

figure()
hold on
plot(rpos, pmaxs, 'x')
plot(r_fit1, y_fit1, ':k', 'LineWidth', 1.5, 'DisplayName', 'Fit: a/r')
xscale('log')
% yscale('log')

%% Save as mat (if desired)
% shower1e9 = struct();
% shower1e9.rpos = rpos;
% shower1e9.p = pmaxs;
% shower1e9.r_fit = r_fit1;
% shower1e9.y_fit = y_fit1;
% save('shower1e9.mat', 'shower1e9');

%% Plot sphere and cylinder with linear and log y-scale
load("cylinder.mat");
load("sphere.mat");
load("shower1e9.mat");  % Load the third dataset (pulse)

p_max_sphere = max(sphere.p);
p_max_cylinder = max(cylinder.p);
p_max_shwer = max(shower1e9.p);  % Normalize third trace

figure;

% === Subplot 1: Linear y-axis ===
subplot(1,2,1);
hold on;

plot(cylinder.rpos, cylinder.p/p_max_cylinder, 'bx', 'DisplayName', 'Cylindrical source');
plot(sphere.rpos, sphere.p/p_max_sphere, 'rx', 'DisplayName', 'Spherical source');
plot(shower1e9.rpos, shower1e9.p/p_max_shwer, 'ko', 'DisplayName', 'Shower $10^9$ GeV');

plot(cylinder.r_fit, cylinder.y_fit/p_max_cylinder, 'b--', ...
    'DisplayName', '$\mathrm{Cylinder\ fit\ } \propto r^{-1/2}$');
plot(sphere.r_fit, sphere.y_fit/p_max_sphere, 'r--', ...
    'DisplayName', '$\mathrm{Sphere\ fit\ } \propto r^{-1}$');

xlabel('$r$ [m]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Peak pressure [a.u.]', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'XScale', 'log', 'FontSize', 14);
title('Linear y-axis', 'Interpreter', 'latex', 'FontSize', 16);
legend('Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 12);
grid on;
axis tight;

% === Subplot 2: Logarithmic y-axis ===
subplot(1,2,2);
hold on;

plot(cylinder.rpos, cylinder.p/p_max_cylinder, 'bx', 'DisplayName', 'Cylindrical source');
plot(sphere.rpos, sphere.p/p_max_sphere, 'rx', 'DisplayName', 'Spherical source');
plot(shower1e9.rpos, shower1e9.p/p_max_shwer, 'ko', 'DisplayName', 'Shower $10^9$ GeV');

plot(cylinder.r_fit, cylinder.y_fit/p_max_cylinder, 'b--', ...
    'DisplayName', '$\mathrm{Cylinder\ fit\ } \propto r^{-1/2}$');
plot(sphere.r_fit, sphere.y_fit/p_max_sphere, 'r--', ...
    'DisplayName', '$\mathrm{Sphere\ fit\ } \propto r^{-1}$');

xlabel('$r$ [m]', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Peak pressure [a.u.]', 'Interpreter', 'latex', 'FontSize', 14);
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 14);
title('Logarithmic y-axis', 'Interpreter', 'latex', 'FontSize', 16);
legend('Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 12);
grid on;
axis tight;

%% ============ THESIS VERSION
figure;
hold on;

set(groot, 'defaultLineLineWidth', 2);
plot(cylinder.rpos, cylinder.p/p_max_cylinder, 'bx', 'DisplayName', 'Cylindrical source', 'MarkerSize', 15);
plot(sphere.rpos, sphere.p/p_max_sphere, 'rx', 'DisplayName', 'Spherical source', 'MarkerSize', 15);
plot(shower1e9.rpos, shower1e9.p/p_max_shwer, 'ko', 'DisplayName', 'Shower $10^9$ GeV', 'MarkerSize', 15);

plot(cylinder.r_fit, cylinder.y_fit/p_max_cylinder, 'b--', ...
    'DisplayName', '$\mathrm{Cylinder\ fit\ } \propto r^{-1/2}$');
plot(sphere.r_fit, sphere.y_fit/p_max_sphere, 'r--', ...
    'DisplayName', '$\mathrm{Sphere\ fit\ } \propto r^{-1}$');

xlabel('$r$ [m]', 'Interpreter', 'latex');
ylabel('Peak pressure [a.u.]', 'Interpreter', 'latex');

set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'FontSize', 30, 'FontName', 'Times New Roman');

% title('Logarithmic y-axis', 'Interpreter', 'latex');
legend('Interpreter', 'latex', 'Location', 'southwest');
grid on;
axis tight;
