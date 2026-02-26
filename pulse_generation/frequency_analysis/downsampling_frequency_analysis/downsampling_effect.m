%% Script to evaluate the compound effects of:
% -> Resampling of a 1 MHz generated signal to 144 kHz
% -> Filtering a signal with the hydrophone transfer function

% Signal gen specs
E0 = 1e10;
r = 500;
z = 10;

addpath('..');
tp_1000 = one_pulse(E0,r,z,0,1e6);
t_1000 = tp_1000(:,1);
p_1000 = tp_1000(:,2);

tp_144  = one_pulse(E0,r,z,0, 144e3);
t_144 = tp_144(:,1); 
p_144 = tp_144(:,2); 

% Filter design (Cecile) -> only for f = 144 kHz
a = [1, -1.07678093, 0.82820418];   % FIR denominator -> frequency dependent!
b = [1, 0, 0];                      % FIR nominator
hydrophoneFIR = dfilt.df2t(b, a);


%% resample & Filter
p_144d = resample(p_1000, 144, 1000);
t_144d = resample(t_1000, 144, 1000);

p_144_filt  = hydrophoneFIR.filter(p_144);
p_144d_filt = hydrophoneFIR.filter(p_144d);

%% plot downsampling effect
ttl = sprintf('Pulse shape for 1 MHz, 144kHz & downsampling at E0 = 1e%d GeV, r = %d m, z = %d m',log10(E0), round(r), round(z));
figure
hold on
title(ttl)
plot(t_1000*1e3, p_1000*1e3, 'DisplayName', 'f = 1 MHz')
plot(t_144*1e3,  p_144*1e3,  'DisplayName', 'f = 144 kHz')
plot(t_144d*1e3, p_144d*1e3, 'DisplayName', 'f = Downsampled 144 kHz')
xlim([-2 2])
set(gca, 'FontSize', 14); % set tick font size
legend
xlabel("t [ms]")
ylabel("p [mPa]")
grid on
hold off

%% plot filtering effect
ttl = sprintf('Pulse shapes at E0 = 1e%d GeV, r = %d m, z = %d m',log10(E0), round(r), round(z));
figure
hold on
title(ttl)  
plot(t_1000*1e3, p_1000*1e3,      'DisplayName', 'f = 1 MHz')
plot(t_144*1e3,  p_144_filt*1e3,  'DisplayName', 'f = 144 kHz, filtered')
plot(t_144d*1e3, p_144d_filt*1e3, 'DisplayName', 'f = Downsampled 144 kHz & filtered')
xlim([-2 2])
set(gca, 'FontSize', 12); % set tick font size
legend
xlabel("t [ms]")
ylabel("p [mPa]")
grid on
hold off