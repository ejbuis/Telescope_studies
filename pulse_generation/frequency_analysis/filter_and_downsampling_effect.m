% script to evaluate the effect of downsampling from any generation
% frequency to 144 kHz, in combination with the effect of the hydrophone
% respone

clear all; close all;
addpath('..');
addpath('../auxiliary');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = 100e6;  % Generation frequency
n_iter = 50; % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fixed parameters
Eo      = 1e10;      
rpos    = 50;       
zpos    = 10;       
nmc     = 1e7;      
c       = 1500;     
atten   = 1;        

% Filter params
G_dB        = 5.486;
G_res_dB    = 28.543;
f_res       = 16200;
fw_res      = 750;

% Sampling frequencies (excluding 1 MHz baseline)
fds = [500e3, 250e3, 144e3]; 
% fds = [144e3, 100e3, 50e3, 32e3]; 
Nfs = length(fds);

% Preallocate retention storage
retention_filt = zeros(n_iter, Nfs);
retention_down = zeros(n_iter, Nfs);  % <-- for downsampled signals

%% Generate traces
fprintf("Generating %d traces for fs = %d MHz \n", n_iter, fs/1e6);
for iter = 1:n_iter
    fprintf("%d ", iter);

    % Generate shower
    Do  = [rpos 0 zpos];        
    rsc = [.5:9.5 15:10:105];      
    zsc = 10:20:2000;              

    tsmc = ShowerParm(rsc, zsc, Eo, 'CORSIKA');
    tsmc = tsmc * diag(kron([1 10], ones(1, 10)));
    pointsc = MCGEn(tsmc, [0 zsc+10], [0 rsc+[0.5*ones(1, 10) 5*ones(1, 10)]],nmc);

    [x,y,z] = pol2cart(rand(nmc, 1)*2*pi, pointsc(:, 2), pointsc(:, 1));
    points = [x y z] * 1e-2; % [m]

    % Generate pulse
    wind = 1024*128;       
    t_axis = (-wind:wind-1)/fs; 

    [~, p, pw, ~] = kernelfrmod(points, Do, log10(Eo), atten, 1, fs, wind);
    fs_orig = fs;

    [b1MHz, a1MHz] = hydrophone(G_dB, G_res_dB, f_res, fw_res, fs_orig);
    p_1MHz_filt = filter(b1MHz, a1MHz, p);

    % Store for example plot (only first iteration)
    if iter == 1
        t_axis_ex = t_axis;
        p_ex = p;
        p_1MHz_filt_ex = p_1MHz_filt;
        p_downsamp_all = cell(1, Nfs);
        t_downsamp_all = cell(1, Nfs);
        p_filt_all = cell(1, Nfs);
    end

    % Downsample + Filter
    for i = 1:Nfs
        f_i = fds(i);
        [b, a] = hydrophone(G_dB, G_res_dB, f_res, fw_res, f_i);

        % Downsample
        p_downsamp = resample(p, f_i, fs_orig);
        t_downsamp = resample(t_axis, f_i, fs_orig);

        % Filter
        p_filt = filter(b, a, p_downsamp);

        % Retention (unfiltered and filtered)
        retention_down(iter, i) = max(p_downsamp) / max(p);              % <-- New
        retention_filt(iter, i) = max(p_filt) / max(p_1MHz_filt);        % Already present

        % Save example trace (first iteration)
        if iter == 1
            p_downsamp_all{i} = p_downsamp;
            t_downsamp_all{i} = t_downsamp;
            p_filt_all{i} = p_filt;
        end
    end
end

%% Compute and display statistics
fprintf('\n===== Amplitude Retention Statistics =====\n');
for i = 1:Nfs
    f_kHz = fds(i)/1e3;

    mean_down = mean(retention_down(:,i));
    std_down  = std(retention_down(:,i));
    mean_filt = mean(retention_filt(:,i));
    std_filt  = std(retention_filt(:,i));

    fprintf('%.0f kHz:\n', f_kHz);
    fprintf('  Downsampled:              Mean = %.3f, Std = %.3f\n', mean_down, std_down);
    fprintf('  Downsampled & Filtered:   Mean = %.3f, Std = %.3f\n', mean_filt, std_filt);
end

%% Plot Example Results (from first iteration)
ttlstring = strcat(round(num2str(fs/1e6)), ' MHz');
colors = lines(Nfs);
figure;
subplot(1,2,1); hold on;
plot(t_axis_ex*1e3, p_ex, 'k-', 'LineWidth', 1.2, 'DisplayName', strcat(round(num2str(fs/1e6)), ' MHz Raw'));
for i = 1:Nfs
    plot(t_downsamp_all{i}*1e3, p_downsamp_all{i}, '-', ...
        'Color', colors(i,:), 'LineWidth', 1.2, ...
        'DisplayName', sprintf('%.0f kHz', fds(i)/1e3));
end
title('Downsampled Pressure Traces (No Filtering)');
xlabel('Time (ms)');
ylabel('Pressure amplitude');
xlim([-2e-1 0]); % ms!
legend('show'); grid on;
fontsize(16, 'points');

subplot(1,2,2); hold on;
plot(t_axis_ex*1e3, p_1MHz_filt_ex, 'k--', 'LineWidth', 1.2, 'DisplayName', strcat(round(num2str(fs/1e6)), ' MHz Filtered'));
for i = 1:Nfs
    plot(t_downsamp_all{i}*1e3, p_filt_all{i}, '--', ...
        'Color', colors(i,:), 'LineWidth', 1.2, ...
        'DisplayName', sprintf('%.0f kHz filtered', fds(i)/1e3));
end
title('Downsampled + Filtered Pressure Traces');
xlabel('Time (ms)');
ylabel('Pressure amplitude');
xlim([-2e-1 0]); % ms!
legend('show'); grid on;
fontsize(16, 'points');
