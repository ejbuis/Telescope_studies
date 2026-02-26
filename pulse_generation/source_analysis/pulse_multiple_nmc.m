% Script to check the effect of the number of monte carlo points on the
% quality of the resulting pulse

addpath('..')
clear all
close all
clc
nmcs = [1e6, 5e6, 1e7, 5e7, 1e8]; % Number of MC points
% nmcs = [1e5,2e5,1e6];

%% Main code
zpos = 20;% m
rpos = 60;% m
Do = [rpos, 0, zpos];
rsc=[.5:9.5 15:10:105]; %radial bin centres (cm)
zsc=10:20:2000;   %Longitudinal Bin Centres (cm)
Eo=1e11; %Primary Energy
atten=1; % Learned's attenuation
fs=1e6; %sampling frequency

%t_axis=(-512:511)/fs; %time axis for plot (default 1024 points)
t_axis=(-1024:1023)/fs; % if you redefine it her, also in other scripts...

amp_dat = zeros(length(nmcs), length(t_axis));

%% Generate energy profile
tsmc = ShowerParm(rsc, zsc, Eo, 'CORSIKA');

%as the 10-100cm bins are 10x wider need to scale by a factor of 10 
tsmc = tsmc * diag(kron([1 10], ones(1,10)));

%% Loop over the different # of Monte-Carlo points
for i=1:length(nmcs)
    nmc = nmcs(i);
    %generate MC points. Note bin EDGES need to be provided 
    pointsc = MCGEn(tsmc, [0 zsc + 10], [0 rsc + [0.5 * ones(1, 10) 5 * ones(1, 10)]], nmc);
    
    %Convert to cartesian
    [x, y, z] = pol2cart(rand(nmc, 1) * 2 * pi, pointsc(:,2), pointsc(:,1));
    % Convert fom cm to m 
    points = [x y z] * 1e-2;
    p = kernelfr2(points, Do, log10(Eo), atten, 10);
    amp_dat(i,:) = p;
end


%%
figure
hold on
colors = lines(length(nmcs)); % Use distinct colors
for i=1:length(nmcs)
    plot(t_axis*1e6, amp_dat(i,:)*1e3, 'DisplayName', ...
        strcat('nmc = ', num2str(nmcs(i)/1e6), ' \times 1e6'), ...
        'LineWidth', 1.5, 'Color', colors(i,:)); % Thicker lines
end
legend('Location', 'northeast', 'FontSize', 12) % Adjust legend position and size
xlabel('Time (\mus)', 'FontSize', 14);
ylabel('Pressure (mPa)', 'FontSize', 14);
grid on
xlim([min(t_axis*1e6), max(t_axis*1e6)]) % Set x-axis limits
ylim([-10, 15]) % Adjust y-axis if needed
set(gca, 'FontSize', 12) % Increase axis font size
title(strcat('E0 = 10^', '{',num2str(log10(Eo)),'}' , ' GeV,  r = ', num2str(rpos), 'm, z = ', num2str(zpos), ' m'))
hold off

%%
numPlots = length(nmcs); 
numRows = ceil(sqrt(numPlots)); 
numCols = ceil(numPlots / numRows); 

figure
title(strcat('E0 = 10^', '{',num2str(log10(Eo)),'}' , ' GeV,  r = ', num2str(rpos), 'm, z = ', num2str(zpos), ' m'))
hold on
for i = 1:numPlots
    subplot(numRows, numCols, i)
    plot(t_axis*1e6, amp_dat(i,:)*1e3, 'LineWidth', 1.5) 
    title(strcat('nmc = ', num2str(nmcs(i)/1e6), ' \times 1e6'))
    xlabel('Time (\mus)', 'FontSize', 12)
    ylabel('Pressure (mPa)', 'FontSize', 12)
    grid on
end

hold off
% set(gcf, 'Position', [100, 100, 1000, 800]) % Resize figure for readability
