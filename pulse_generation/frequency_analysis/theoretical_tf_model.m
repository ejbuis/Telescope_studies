%% Script to generate theoretical tf model:
% STEP 0: Load data
clear all; close all; clc
addpath('..');

% Response at 90 degrees for a double membrane hydrophone
load('complex_hydrophone_response.mat'); % -> loads fs and H_complex
H_dB = 20 * log10(abs(H_complex));      % attenuation measurement
H_phase = unwrap(angle(H_complex));     % phase measurement
fs_filt = 1e6;                          % sampling frequenct of traces

% Parameters (in dB) - 90 degrees
G_dB = 5;                % Baseline gain in dB (e.g., 0 dB = unity gain)
G_res_dB = 28.54;        % Peak gain at resonance in dB
f_res = 16210;           % Resonance frequency (Hz)
fw_res = 750;            % Full width at half maximum (Hz)

[b,a] = hydrophone(G_dB, G_res_dB, f_res, fw_res, fs_filt);
H = tf(b, a, 1/fs_filt);
% figure;
% bode(H);

N_fft = 8192;
[H_filt, f_filt] = freqz(b, a, N_fft, fs_filt);  % f_filt is in Hz

% Convert to dB and unwrap phase
H_filt_dB = 20 * log10(abs(H_filt));
H_filt_phase = unwrap(angle(H_filt));

% Plot magnitude
figure;
subplot(2,1,1);
hold on;
semilogx(fs, H_dB, 'b', 'DisplayName', 'Measured', 'LineWidth', 1.5);
semilogx(f_filt, H_filt_dB, 'r--','DisplayName', 'Theoretical model', 'LineWidth', 1.5);
xscale('log');
xlim([min(fs) max(fs)]);
ylabel('Magnitude (dB)');
grid on;
legend;
title('FIR Fit (Magnitude and Phase)');
fontsize(16, 'points')

% Plot phase
subplot(2,1,2);
hold on;
xscale('log');
xlim([min(fs) max(fs)]);
semilogx(fs, rad2deg(H_phase), 'b', 'DisplayName', 'Measured', 'LineWidth', 1.5);
semilogx(f_filt, rad2deg(H_filt_phase), 'r--', 'DisplayName', 'Theoretical model', 'LineWidth', 1.5);
ylabel('Phase (degrees)');
fontsize(16, 'points')
xlabel('Frequency (Hz)');
grid on;
legend;

%% angular sensitivity, with combined subplot
% Use the spectra generated at the various angles

angls = [0,30,60,90,120,150,180];
G_dBs       = zeros(size(angls));
G_res_dBs   = zeros(size(angls));
f_ress      = zeros(size(angls));
fw_ress     = zeros(size(angls));

% Limit the search of peaks
f_peak_min = 14e3;
f_peak_max = 18e3;
i_peak_min = find(fs > f_peak_min, 1, 'first');
i_peak_max = find(fs < f_peak_max, 1, 'last');

% Limit the search of steady-state
f_ss_min = 8e3;
f_ss_max = 10e3;
i_ss_min = find(fs > f_ss_min, 1, 'first');
i_ss_max = find(fs < f_ss_max, 1, 'last');

% Create one figure with two subplots
figure;
subplot(2,1,1); % Magnitude plot
hold on;
grid on;
ylabel('Magnitude (dB)');
title('FIR Fit Magnitude at Various Angles');
fontsize(16, 'points')
xlim([min(fs) max(fs)]);

subplot(2,1,2); % Phase plot
hold on;
grid on;
ylabel('Phase (degrees)');
xlabel('Frequency (Hz)');
title('FIR Fit Phase at Various Angles');
fontsize(16, 'points')
xlim([min(fs) max(fs)]);

% Loop over all angles and add traces to both subplots
for i = 1:length(angls)
    fname = strcat('./Transfer_functions/response',num2str(angls(i)), 'deg.mat');
    load(fname);

    % Original data for comparison
    H_dB = 20 * log10(abs(H_complex));
    H_phase = unwrap(angle(H_complex));

    % Get parameters
    [G_res_dB_i, indmax] = max(H_dB(i_peak_min:i_peak_max));
    G_dB_i = mean(H_dB(i_ss_min:i_ss_max));

    f_res_i = fs(indmax+i_peak_min-1);
    fw_res_i = fw_res;
    
    G_dBs(i)    = G_dB_i;   G_res_dBs(i)    = G_res_dB_i;
    f_ress(i)   = f_res_i;  fw_ress(i)      = fw_res_i;
    [b_i, a_i] = hydrophone(G_dB_i, G_res_dB_i, f_res_i, fw_res_i, fs_filt);

    [H_filt_i, f_filt_i] = freqz(b_i, a_i, N_fft, fs_filt);  % f_filt is in Hz
    
    H_dB_i = 20 * log10(abs(H_filt_i));
    H_phase_i = unwrap(angle(H_filt_i));

    % Plot magnitude
    subplot(2,1,1);
    h_mag = semilogx(fs, H_dB, 'DisplayName', sprintf('%d°', angls(i)));
    semilogx(f_filt_i, H_dB_i, '--', 'Color', h_mag.Color, 'HandleVisibility', 'off');

    % Plot phase
    subplot(2,1,2);
    h_phase = semilogx(fs, rad2deg(H_phase), 'DisplayName', sprintf('%d°', angls(i)));
    semilogx(f_filt_i, rad2deg(H_phase_i), '--', 'Color', h_phase.Color, 'HandleVisibility', 'off');
end

% Add legends after the loop
subplot(2,1,1);
xline(f_ss_min, 'r-', 'HandleVisibility', 'off');
xline(f_peak_max, 'r-', 'HandleVisibility', 'off');
legend;

subplot(2,1,2);
xline(f_ss_min, 'r-', 'HandleVisibility', 'off');
xline(f_peak_max, 'r-', 'HandleVisibility', 'off');
legend;

%% Thesis version of the plot above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create one figure with two subplots

angls_plot = [30,90,120];
figure;
hold on;
grid on;
ylabel('Response magnitude [dB re \muPa/\muPa]');
xlabel('Frequency [kHz]')
% title('FIR Fit Magnitude at Various Angles');
xlim([min(fs)/1000 max(fs)/1000]);
set(gca, 'FontSize', 30, 'FontName', 'Times New Roman');
set(groot, 'defaultLineLineWidth', 3);
set(gca, 'GridLineWidth',1.5)

% Loop over all angles and add traces to both subplots
for i = 1:length(angls_plot)
    fname = strcat('./Transfer_functions/response',num2str(angls_plot(i)), 'deg.mat');
    load(fname);

    % Original data for comparison
    H_dB = 20 * log10(abs(H_complex));
    H_phase = unwrap(angle(H_complex));

    % Get parameters
    [G_res_dB_i, indmax] = max(H_dB(i_peak_min:i_peak_max));
    G_dB_i = mean(H_dB(i_ss_min:i_ss_max));

    f_res_i = fs(indmax+i_peak_min-1);
    fw_res_i = fw_res;
    
    G_dBs(i)    = G_dB_i;   G_res_dBs(i)    = G_res_dB_i;
    f_ress(i)   = f_res_i;  fw_ress(i)      = fw_res_i;
    [b_i, a_i] = hydrophone(G_dB_i, G_res_dB_i, f_res_i, fw_res_i, fs_filt);

    [H_filt_i, f_filt_i] = freqz(b_i, a_i, N_fft, fs_filt);  % f_filt is in Hz
    
    H_dB_i = 20 * log10(abs(H_filt_i));
    H_phase_i = unwrap(angle(H_filt_i));

    h_mag = semilogx(fs / 1000, H_dB, 'DisplayName', sprintf('θ = %d°', angls_plot(i)));
    semilogx(f_filt_i / 1000, H_dB_i, '--', 'Color', h_mag.Color, 'HandleVisibility', 'off');

end


xline(f_ss_min / 1000, 'k-', 'DisplayName', 'fit range', 'LineWidth', 2);
xline(f_peak_max / 1000, 'k-', 'HandleVisibility', 'off', 'LineWidth', 2);
legend('Interpreter', 'latex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot resonance peaks as function of angles
G_res_lin = 10.^(G_res_dBs./20);

figure;
hold on;
plot(angls, G_res_lin, 'x-', 'LineWidth', 1.5, 'MarkerSize', 10);

% Axes labels
xlabel('Angle [°]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Gain [-]', 'FontSize', 12, 'FontWeight', 'bold');

% Title
title('Linear gain at resonance vs. angle', 'FontSize', 14, 'FontWeight', 'bold');

% Grid and box
grid on;
box on;

% X axis ticks in degrees
xticks(0:30:180); % adjust range if needed
xlim([min(angls) max(angls)]);

% Y axis formatting
ylim([0 max(G_res_lin)*1.1]);

% Improve appearance
set(gca, 'FontSize', 12, 'LineWidth', 1);
