%% Generation frequency script enhanced with statistics on multiple runs

clear all; close all;
addpath('../..');
addpath('../../auxiliary');

% Parameters
n_iter = 50;
Eo     = 1e10;
rpos   = 50;
zpos   = 10;
nmc    = 1e7;
atten  = 1;
wind   = 1024*128;
fontsize = 14;

fss = logspace(3, 8, 50);   % 1 kHz to 100 MHz
fss = [500e3, 750e3, 1e6, 5e6, 10e6, 25e6, 50e6, 75e6, 100e6];
Nfreqs = length(fss);

% Preallocate
pmaxoffsets_all = zeros(n_iter, Nfreqs);
pminoffsets_all = zeros(n_iter, Nfreqs);
opss_all        = zeros(n_iter, Nfreqs);
comptimes_all   = zeros(n_iter, Nfreqs);

% Storage for iteration 1
trace_storage = struct();

% Progress bar
h = waitbar(0, 'Running iterations...');

for iter = 1:n_iter
    waitbar(iter/n_iter, h, sprintf('Running iteration %d of %d', iter, n_iter));

    % === Monte Carlo Event Generation ===
    rsc  = [.5:9.5 15:10:105];       
    zsc  = 10:20:2000;
    tsmc = ShowerParm(rsc, zsc, Eo, 'CORSIKA');
    tsmc = tsmc * diag(kron([1 10], ones(1, 10)));
    pointsc = MCGEn(tsmc, [0 zsc+10], ...
              [0 rsc+[0.5*ones(1, 10) 5*ones(1, 10)]], nmc);
    [x, y, z] = pol2cart(rand(nmc, 1)*2*pi, pointsc(:,2), pointsc(:,1));
    points = [x y z] * 1e-2;  % m
    Do = [rpos 0 zpos];

    pmaxs = zeros(1, Nfreqs);
    pmins = zeros(1, Nfreqs);
    opss  = zeros(1, Nfreqs);
    comptimes = zeros(1, Nfreqs);

    for i = 1:Nfreqs
        fs = fss(i);
        t_axis = (-wind:wind-1)/fs;

        tic;
        [~, p, pw, ~] = kernelfrmod(points, Do, log10(Eo), atten, 1, fs, wind);
        comptimes(i) = toc;

        pmaxs(i) = max(p);
        pmins(i) = min(p);
        opss(i)  = trapz(abs(p.^2)) * (t_axis(2) - t_axis(1));

        if iter == 1
            trace_storage.ts(i,:) = t_axis;
            trace_storage.ps(i,:) = p;
            f = (0:(length(pw)/2-1)) * (fs/length(pw));
            trace_storage.pws(i,:) = abs(pw(1:length(f))).^2;
            trace_storage.f_pwss(i,:) = f;
        end
    end

    % Normalize w.r.t. highest fs
    pmaxoffsets_all(iter,:) = pmaxs / pmaxs(end);
    pminoffsets_all(iter,:) = pmins / pmins(end);
    opss_all(iter,:)        = opss;
    comptimes_all(iter,:)   = comptimes;
end
close(h);

%% === Stats ===
mean_max = mean(pmaxoffsets_all, 1);
std_max  = std(pmaxoffsets_all, 0, 1);
mean_min = mean(pminoffsets_all, 1);
std_min  = std(pminoffsets_all, 0, 1);
fss_kHz  = fss * 1e-3;

%% === Plot 1: Pressure Traces (from iteration 1) ===
figure;
subplot(1,2,1); hold on
for i = 1:Nfreqs
    plot(trace_storage.ts(i,:), trace_storage.ps(i,:) * 1e3, ...
         'DisplayName', sprintf('f = %.0f kHz', fss_kHz(i)));
end
xlabel('Time (s)', 'FontSize', fontsize);
ylabel('Pressure (mPa)', 'FontSize', fontsize);
xlim([-2e-4 2e-4])
title(sprintf('Pressure traces (1 example)'), 'FontSize', fontsize);
legend; grid on; set(gca, 'FontSize', fontsize); hold off

subplot(1,2,2); hold on
for i = 1:Nfreqs
    plot(trace_storage.f_pwss(i,:), trace_storage.pws(i,:), ...
         'DisplayName', sprintf('f = %.0f kHz', fss_kHz(i)));
end
xlabel('Frequency (Hz)', 'FontSize', fontsize);
ylabel('Power', 'FontSize', fontsize);
title('Power Spectrum', 'FontSize', fontsize);
legend; grid on; set(gca, 'FontSize', fontsize); hold off

%% === Plot 2: Retention + Errorbars ===
ttl = sprintf('Amplitude retention for E0 = 1e%d GeV, r = %d m, z = %d m, nmc = 1e%.2f', ...
               log10(Eo), round(rpos), round(zpos), log10(nmc));
figure; hold on
title(ttl, 'FontSize', fontsize)
errorbar(fss_kHz, mean_max, std_max, 'x-', 'LineWidth', 1.5)
errorbar(fss_kHz, mean_min, std_min, 'o-', 'LineWidth', 1.5)
yline(1, '--', '100% retention')
xscale('log');
ylim([0.4 1.1]);
xlabel("Signal gen frequency [kHz]"); ylabel("p_{f} / p_{f=100MHz}")
legend('max pressure', 'min pressure')
set(gca, 'FontSize', fontsize); grid on; hold off

%% === Plot 3: Integral (1 example, could also use mean ± std) ===
figure; hold on
title('Integral value of |p|', 'FontSize', fontsize)
plot(fss_kHz, mean(opss_all,1), 'x-', 'MarkerSize', 8)
xscale('log')
xlabel("Signal gen frequency [kHz]"); ylabel("\Sigma |p^2| dt")
set(gca, 'FontSize', fontsize); grid on; hold off

%% === Plot 4: Computation time ===
figure; hold on
title('Computation time vs frequency', 'FontSize', fontsize)
plot(fss_kHz, mean(comptimes_all,1), 'x-', 'MarkerSize', 8)
xscale('log')
xlabel("Signal gen frequency [kHz]"); ylabel("Computation time [s]")
set(gca, 'FontSize', fontsize); grid on; hold off
