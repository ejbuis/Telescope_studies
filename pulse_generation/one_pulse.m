function tp = one_pulse(Eo, rpos, zpos, plot_bool, fs, name, bin)             
% Script to generate 1 pulse & save to a file from r,z specification
% Eo            Primary Energy        [GeV],  default 1e9 GeV
% rpos          Observer r position   [m],    default 1000 m
% zpos          Observer z position   [m],    default 6 m
% plot_bool     Whether to plot fit   [-],    default 1
% fs            Sampling frequency    [Hz],   default 1e6 Hz
% name          Name to save file     [-],    if undefined, not saved
% bin           If set to 1 save as binary
if nargin < 1 Eo = 1e9;         end
if nargin < 2 rpos = 1e3;       end
if nargin < 3 zpos = 6;         end
if nargin < 4 plot_bool = 1;    end
if nargin < 5 fs = 1e6;         end
if nargin < 6 save_bool = 0;    else save_bool = 1; end
if nargin < 7 bin = 0;          end

tic
%% Main code
rsc=[.5:9.5 15:10:105]; %radial bin centres (cm)
zsc=10:20:2000;   
Do=[rpos 0 zpos]; % Position of observer         

%t_axis=(-512:511)/fs; %time axis for plot (default 1024 points)
wind = 4096*4;

t_axis=(-wind:(wind-1))/fs; % if you redefine it her, also in other scripts...
atten=5; % Learned's attenuation = 1, No attenuation = 5
% disp("USING NO ATTENUATION");
nmc=2e6; % Number of MC points 
tsmc=ShowerParm(rsc,zsc,Eo,'CORSIKA');

%as the 10-100cm bins are 10x wider need to scale by a factor of 10 
tsmc=tsmc*diag(kron([1 10],ones(1,10)));

%generate MC points. Note bin EDGES need to be provided 
pointsc=MCGEn(tsmc,[0 zsc+10],[0 rsc+[0.5*ones(1,10) 5*ones(1,10)]],nmc);

%Convert to cartesian
[x,y,z]=pol2cart(rand(nmc,1)*2*pi,pointsc(:,2),pointsc(:,1));

% Convert fom cm to m 
points=[x y z]*1e-2;

% [p, pw, E]=kernelfr2(points,Do,log10(Eo),atten,10,fs);
disp("USING KERNELFRMOD");
[t_ret, p, pw, E]=kernelfrmod(points,Do,log10(Eo),atten,10,fs,wind);
            

% save the data 
%filename = args{2}
tp = [t_axis' p]; % check this weird construction of the transpose columns
if save_bool
    if bin
        filename = strcat(name, "_", num2str(zpos), "_", num2str(rpos));
        fid = fopen(filename, 'w');       % open file for writing (binary)
        fwrite(fid, p, 'double');         % write p as double-precision binary
        fclose(fid);
    else
        filename = strcat(name, "_", num2str(zpos), "_", num2str(rpos)); %, ".dat");
        save(filename, 'tp', '-ascii'); % close the file
    end
end

%% Power spectrum
N = length(pw);               % Length of FFT
f = (0:N/2-1)*(fs/N);         % Frequency axis (only positive frequencies)
power_spectrum = abs(pw(1:N/2)).^2;  % Compute power spectrum (one-sided)
toc

%% Plotting
if plot_bool
    figure;
    hold on
    plot(t_axis*1e6,p*1000);
    % plot(t_axis*1e6,pfilt*1000);
    xlabel('Time (\mus)');
    ylabel('Pressure (mPa)');
    % xlim([-100 100])
    % ylim([-12 12])
    grid on
    hold off

    fontsize(16, 'points');
    
    % Plot
    figure;
    plot(f*1e-3, power_spectrum/max(power_spectrum), 'LineWidth', 2);
    xlabel('Frequency (kHz)');
    xlim([0,72]);
    xline(20, '--k', 'LineWidth', 2);
    xline(16, '--r', {'hydrophone','resonance'}, 'LineWidth', 2, 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');
    ylabel('Power [a.u.]');
    % yscale('log');
    title('Power Spectrum');
    grid on;
    fontsize(12, 'points');
end

end
