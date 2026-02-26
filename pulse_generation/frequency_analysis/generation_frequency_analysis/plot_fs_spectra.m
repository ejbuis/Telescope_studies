% script used to generate the pulses_and_spectra.jp plot

addpath('..');
addpath('../..');
Eo = 1e11;
wind = 8192;
rpos = 1000;
zpos = 6;

%% Set params
rsc     = [.5:9.5 15:10:105]; %radial bin centres (cm)
zsc     = 10:20:2000;   
Do      = [rpos 0 zpos]; % Position of observer      
atten   = 1; % Learned's attenuation = 1, No attenuation = 5
nmc     = 2e6; % Number of MC points 
tsmc    = ShowerParm(rsc,zsc,Eo,'CORSIKA');

%as the 10-100cm bins are 10x wider need to scale by a factor of 10 
tsmc    = tsmc*diag(kron([1 10],ones(1,10)));

%generate MC points. Note bin EDGES need to be provided 
pointsc = MCGEn(tsmc,[0 zsc+10],[0 rsc+[0.5*ones(1,10) 5*ones(1,10)]],nmc);

%Convert to cartesian
[x,y,z] = pol2cart(rand(nmc,1)*2*pi,pointsc(:,2),pointsc(:,1));

% Convert fom cm to m 
points  = [x y z]*1e-2;

%%
[t_1MHz, p_1MHz,pw_1MHz,~]          = kernelfrmod(points,Do,log10(Eo),atten,10,1e6,wind);

[t_144kHz, p_144kHz,pw_144kHz,~]    = kernelfrmod(points,Do,log10(Eo),atten,10,144e3,wind);

p_144kHz_resamp = resample(p_1MHz, 144, 1000, 2);
scaling = 1; % max(p_1MHz) / max(p_144kHz_resamp);
p_144kHz_resamp = p_144kHz_resamp * scaling; % rescale
t_144kHz_resamp = resample((t_1MHz - min(t_1MHz)), 144, 1000, 2) + min(t_1MHz);

t_1MHz = t_1MHz - t_1MHz(wind);
t_144kHz = t_144kHz - t_144kHz(wind);
t_144kHz_resamp = t_144kHz_resamp - t_144kHz_resamp(length(t_144kHz_resamp)/2) - (t_144kHz_resamp(2)-t_144kHz_resamp(1))/2;
%% Power spectrum
pp_1MHz             = fft(p_1MHz*1e6);
pp_144kHz           = fft(p_144kHz*1e6);
pp_144kHz_resamp    = fft(p_144kHz_resamp*1e6);

% Length of FFT
N1MHz               = length(pp_1MHz);               
N144kHz             = length(pp_144kHz);         
N144kHz_resamp      = length(pp_144kHz_resamp);      

% Frequency axis (only positive frequencies)
f1MHz               = (0:N1MHz/2-1)*(1e6/N1MHz);         
f144kHz             = (0:N144kHz/2-1)*(144e3/N144kHz);     
f144kHz_downsamp    = (0:N144kHz_resamp/2-1)*(144e3/N144kHz_resamp);     

% Compute power spectrum (one-sided)
P_1MHz              = abs(pp_1MHz(1:N1MHz/2) / N1MHz).^2;  
P_144kHz            = abs(pp_144kHz(1:N144kHz/2) / N144kHz).^2;  
P_144kHz_resamp     = abs(pp_144kHz_resamp(1:N144kHz_resamp/2) / N144kHz_resamp).^2;  

%% plot
figure;

% --- Left panel: Pulses ---
subplot(1,2,1); hold on;

plot(t_1MHz*1e3, p_1MHz, 'LineWidth', 1.5);  % 1 MHz pulse
plot(t_144kHz*1e3, p_144kHz, 'LineWidth', 1.5);  % 144 kHz pulse
plot(t_144kHz_resamp*1e3, p_144kHz_resamp, 'LineWidth', 1.5);  % Resampled pulse

xlabel('Time [ms]');
xlim([-0.1 0.1])
ylabel('Amplitude [Pa]');
title('Pulses in Time Domain');
legend('1 MHz', '144 kHz', 'Downsampled');
set(gca, 'FontSize', 14); % set tick font size
grid on;

% --- Right panel: Power spectra ---
subplot(1,2,2); hold on;

plot(f1MHz/1e3, 10*log10(P_1MHz), 'LineWidth', 1.5); % 10*log10(P_1MHz) P_1MHz
plot(f144kHz/1e3, 10*log10(P_144kHz), 'LineWidth', 1.5); % 10*log10(P_144kHz) P_144kHz
plot(f144kHz_downsamp/1e3, 10*log10(P_144kHz_resamp), 'LineWidth', 1.5); % 10*log10(P_144kHz_resamp) P_144kHz_resamp

xlabel('Frequency [kHz]');
ylabel('Power [dB re \mu Pa]');
title('Power Spectra');
xscale('log')
xlim([0 500]);
xline(1, '--k',label='1 kHz');
xline(25, '--k',label='25 kHz');
xline(72, 'k',label='Nyquist');
set(gca, 'FontSize', 14); % set tick font size
legend('1 MHz', '144 kHz', 'Downsampled');
grid on;




