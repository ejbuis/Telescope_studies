% Script to analyse the effect of a differing sampling frequency on the
% pulse shape

% clear all; close all;
addpath('../..');
addpath('../../auxiliary');

rng(10); % FIX RANDOM SEED!

Eo      = 1e11;         % Primary Energy        [GeV],  only for Corsika
rpos    = 1000;          % Observer r position   [m],    default 1000 m
zpos    = 6;           % Observer z position   [m],    default 6 m
nmc     = 1e6;          % Number of MC points   [-],    default 1e6
c       = 1500;         % speed of sound        [m/s]
atten   = 1;            % Learned's attenuation = 1, No attenuation = 5

tic

%% ============== Energy distribution & MC generation ==============
Do  = [rpos 0 zpos];        % Position of observer

fprintf("==== Generating from CORSIKA Parametrisation ==== \n")

% => Binning specifics
rsc = [.5:9.5 15:10:105];       % Radial bin centres (cm)
zsc = 10:20:2000;               % Longitudinal Bin Centres (cm)

tsmc = ShowerParm(rsc, zsc, Eo, 'CORSIKA');

% as the 10-100cm bins are 10x wider need to scale by a factor of 10 
tsmc = tsmc * diag(kron([1 10], ones(1, 10)));
pointsc = MCGEn(tsmc, [0 zsc+10], [0 rsc+[0.5*ones(1, 10) 5*ones(1, 10)]],nmc);

% Convert to cartesian
[x,y,z] = pol2cart(rand(nmc, 1)*2*pi, pointsc(:, 2), pointsc(:, 1));

% Convert fom cm to m 
points = [x y z] * 1e-2;

%% ============== define frequency axes ============== %% 
% fss = logspace(5,6.301,20);     % range from 100 kHz to 2 MHz
% fss = logspace(5,7,20);       % range from 100 kHz to 10 MHz
fss = logspace(3,8,50);         % range from 1 kHz to 100 MHz
% fss = [144e3, 250e3,500e3,1e6,2e6,5e6];

Nfreqs = length(fss);
opss = zeros(size(fss));        % Array to keep track of the integral value of the pulses
comptimes = zeros(size(fss));   % Array to keep track of the computation times

wind = 1024*128;        % 1/2 size of the time window [-]

ts      = zeros(Nfreqs, wind*2);
ps      = zeros(Nfreqs, wind*2);
pws     = zeros(Nfreqs, wind);
f_pwss  = zeros(Nfreqs, wind);

%% ============== Compute pressure signal ==============
comptimes = zeros(1, Nfreqs); % Preallocate

fprintf("Generating %d signals \n", Nfreqs);
for i = 1:Nfreqs
    fprintf("%d ", i);
    fs = fss(i);
    tic; % Start timing

    t_axis = (-wind:wind-1)/fs; % time axis (centered at 0)
    ts(i,:) = t_axis;

    [~, p, pw, ~] = kernelfrmod(points, Do, log10(Eo), atten, 1, fs, wind);
    ps(i,:) = p;

    opss(i) = trapz(abs(p.^2))*(t_axis(2)-t_axis(1)); % Simple integral

    N = length(pw);                     % Length of FFT
    f = (0:N/2-1) * (fs/N);             % Frequency axis (only positive frequencies)
    power_spectrum = abs(pw(1:N/2)).^2; % Compute power spectrum (one-sided)

    pws(i,:)    = power_spectrum;
    f_pwss(i,:) = f;

    comptimes(i) = toc; % Store time in seconds
end


%% Combine into 1 plot
plot_title = sprintf('Pressure plot for r = %d m, z = %d m, nmc = %.0e', round(rpos), round(zpos), nmc);
fssnice = round(fss/1e3);
fontsize = 14; % desired font size

figure;

% Left subplot: pressure vs time
subplot(1,2,1);
hold on
for i = 1:Nfreqs
    plot(ts(i,:), ps(i,:)*1e3, 'DisplayName', strcat("f = ", num2str(fssnice(i)), " kHz"));
end
xlabel('Time (s)', 'FontSize', fontsize);
ylabel('Pressure (mPa)', 'FontSize', fontsize);
xlim([-2e-4 2e-4])
title(plot_title, 'FontSize', fontsize);
legend;
grid on
set(gca, 'FontSize', fontsize); % set tick font size
hold off

% Right subplot: power spectrum
subplot(1,2,2);
hold on
for i = 1:Nfreqs
    plot(f_pwss(i,:), pws(i,:), 'DisplayName', strcat("f = ", num2str(fssnice(i)), " kHz"));
end
xlabel('Frequency (kHz)', 'FontSize', fontsize);
ylabel('Power', 'FontSize', fontsize);
title('Power Spectrum', 'FontSize', fontsize);
legend;
grid on
set(gca, 'FontSize', fontsize); % set tick font size
hold off

%% Compute and plot max offsets
pmaxs = max(ps,[],2);
pmins = min(ps,[],2);
pmaxoffsets = pmaxs/pmaxs(end);
pminoffsets = pmins/pmins(end);

% Fit exponential retention to data traces.
pfitmax = fit_exp_retention(fss, pmaxoffsets);
pfitmin = fit_exp_retention(fss, pminoffsets);
fssfit = linspace(min(fss),max(fss),100);


ttl = sprintf('Amplitude retention for E0 = 1e%d GeV, r = %d m, z = %d m, nmc = 1e%.2f',log10(Eo), round(rpos), round(zpos), log10(nmc));

figure;
hold on
title(ttl)
plot(fss*1e-3, pmaxoffsets, 'x', 'MarkerSize',10)
plot(fss*1e-3, pminoffsets, 'o', 'MarkerSize',10)
% plot(fssfit*1e-3, pfitmax(fssfit), '--k')
% plot(fssfit*1e-3, pfitmin(fssfit), '--k')
yline(1,'--', '100% retention')
ylim([0.4 1.1])
xscale('log')
xlabel("Signal gen frequency [kHz]")
ylabel("p_{f} / p_{f=100MHz}")
legend('maximum pressure', 'minimum pressure')
set(gca, 'FontSize', fontsize); % set tick font size
grid on
hold off

%% Plot the integral values

ttl2 = sprintf('Integral value of |p|');

figure;
hold on
title(ttl2)
plot(fss*1e-3, opss, 'x', 'MarkerSize',10)
xscale('log')
xlabel("Signal gen frequency [kHz]")
ylabel("\Sigma |p^2|/dt")
set(gca, 'FontSize', fontsize); % set tick font size
grid on
hold off

%% Plot the computation times
ttl3 = sprintf('Computation time as a function of f');

figure;
hold on
title(ttl3)
plot(fss*1e-3, comptimes, 'x', 'MarkerSize',10)
xscale('log')
xlabel("Signal gen frequency [kHz]")
ylabel("Computation time [s]")
set(gca, 'FontSize', fontsize); % set tick font size
grid on
hold off