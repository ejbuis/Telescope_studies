%% This script generates a 144 kHz template bank
% To be used in Jpp implementations

% DO NOT CHANGE:
addpath('..');
addpath('../auxiliary');
c           = 1500;         % [m/s] speed of sound
fs_gen      = 1e6;          % [Hz] signal gen frequency
fs_JPP      = 144e3;        % [Hz] JPP sample frequency
wind        = 4096*2;       % 1/2 size of the time window [-]
Nrot        = 100;          % number of rotations to consider

flpf        = 100;         % low pass filter pass band [Hz]

Norig       = wind * 2;
p_temp      = resample(ones(Norig,1), fs_JPP, fs_gen, 2);
N_resampeld = size(p_temp, 1);
clear p_temp;

%% Extend with more template ideas if desired!
% DO CHANGE:
do_plot = 1;

templates = struct();
templates.template1 = {'near_field_in_pancake', 1e11, 100,  6,  zeros(Norig, 1), zeros(N_resampeld, 1), 0};
templates.template2 = {'mid_field_in_pancake',  1e11, 1000, 6,  zeros(Norig, 1), zeros(N_resampeld, 1), 0};
templates.template3 = {'far_field_in_pancake',  1e11, 2500, 6,  zeros(Norig, 1), zeros(N_resampeld, 1), 0};

templates.template4 = {'far_field_north', 1e11, 2500,  -30,  zeros(Norig, 1), zeros(N_resampeld, 1), 0};
templates.template5 = {'far_field_south', 1e11, 2500, 50,  zeros(Norig, 1), zeros(N_resampeld, 1), 0};

templates.template6 = {'near_field_north',  1e11, 100, -30,  zeros(Norig, 1), zeros(N_resampeld, 1), 1};
templates.template7 = {'near_field_south',  1e11, 100, 50,  zeros(Norig, 1), zeros(N_resampeld, 1), 1};

templates.template8 = {'mid_field_north',  1e11, 1000, -30,  zeros(Norig, 1), zeros(N_resampeld, 1), 0};
templates.template9 = {'mid_field_south',  1e11, 1000, 50,  zeros(Norig, 1), zeros(N_resampeld, 1), 0};

templates.template10 = {'on_axis',  1e11, 1, 6,  zeros(Norig, 1), zeros(N_resampeld, 1), 1};

templateNames = fieldnames(templates);


%% Create shower
rsc = [.5:9.5 15:10:105]; % Radial bin centres (cm)
zsc = 10:20:2000;         % Longitudinal Bin Centres (cm)
atten = 7;                % Learned's attenuation = 1, ACORNE = 7

%% loop over the templates
tic
for i = 1:numel(templateNames)
    
    % unpack template information
    name    = templateNames{i};
    tpl     = templates.(name);
    label   = tpl{1};
    E0      = tpl{2};
    rpos    = tpl{3};
    zpos    = tpl{4};
    dolpf   = tpl{7};

    fprintf("Generating template %i, %s\n", i,  label);

    nmc     = nmc_from_r(rpos);

    Do = [rpos 0 zpos]; % Position of observer

    % NOTE: The following only works for the standard binning as included in
    % the Corsika parameterisation in Showerparm
    tsmc = interpolate_Sm(E0);
    
    % as the 10-100cm bins are 10x wider need to scale by a factor of 10 
    tsmc = tsmc * diag(kron([1 10], ones(1,10)));

    % Determine the monte carlo points, nmc = f(r)
    pointsc = MCGEn(tsmc, [0 zsc + 10], [0 rsc + [0.5 * ones(1, 10) 5 * ones(1, 10)]], nmc);

    % Convert to cartesian
    [x_r, y_r, z_r] = pol2cart(rand(nmc, 1) * 2 * pi, pointsc(:,2), pointsc(:,1));

    % Convert fom cm to m 
    points = [x_r y_r z_r] * 1e-2;

    % Compute the pressure
    [~,p,~,~] = kernelfrmod(points, Do, log10(E0), atten, Nrot, fs_gen, wind);
    pdsamp    = resample(p, fs_JPP, fs_gen, 2);

    if dolpf
        pn = lowpass(p,flpf,fs_JPP);
        figure;
        hold on
        plot(p)
        plot(pn)
        hold off
        p = pn;
    end

    tpl{5} = p;
    tpl{6} = pdsamp;
    templates.(name) = tpl;
end
toc

%% Writing
% === Create unique output folder ===
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
outDir = sprintf('../../template_banks/template_bank_%s', timestamp);
mkdir(outDir);

% === Prepare reference text file ===
listFile = fullfile(outDir, 'template_list.txt');
fid_list = fopen(listFile, 'w');
if fid_list == -1
    error('Could not create template list file.');
end

% Write header with fixed-width columns
fprintf(fid_list, '%-30s %-12s %-12s %-12s\n', 'Name', 'E0 [Hz]', 'rpos [m]', 'zpos [m]');
fprintf(fid_list, '%s\n', repmat('-', 1, 70));

% === Loop through templates ===
templateNames = fieldnames(templates);

for i = 1:numel(templateNames)
    name = templateNames{i};
    tpl = templates.(name);

    % Extract parameters
    label = tpl{1};
    E0    = tpl{2};
    rpos  = tpl{3};
    zpos  = tpl{4};
    data  = tpl{6};

    % Sanitize filename
    safe_label = regexprep(label, '[^\w\-]', '_');
    binFile = fullfile(outDir, sprintf('%s.bin', safe_label));

    % Write binary data
    fid_bin = fopen(binFile, 'w');
    if fid_bin == -1
        error('Could not open %s for writing.', binFile);
    end
    fwrite(fid_bin, data, 'double');
    fclose(fid_bin);

    % Write aligned line in reference list
    fprintf(fid_list, '%-30s %-12.3e %-12.3g %-12.3g\n', label, E0, rpos, zpos);
end

% === Close file and print info ===
fclose(fid_list);

fprintf('Templates saved in folder: %s\n', outDir);
fprintf('Reference file: %s\n', listFile);



%% Plotting
if do_plot
    figure;
    hold on;
    
    for i = 1:numel(templateNames)
        name = templateNames{i};
        tpl = templates.(name);
        
        data = tpl{5};  % The vector you want to plot
        plot(data, 'DisplayName', tpl{1});  % tpl{1} is the descriptive label
    end
    
    hold off;
    legend('show');
    xlabel('Sample index');
    ylabel('Amplitude');
    title('ORIGINAL signals');
    grid on;
    
    figure;
    hold on;
    
    for i = 1:numel(templateNames)
        name = templateNames{i};
        tpl = templates.(name);
        
        data = tpl{6};  % The vector you want to plot
        plot(data, 'DisplayName', tpl{1});  % tpl{1} is the descriptive label
    end
    
    hold off;
    legend('show');
    xlabel('Sample index');
    ylabel('Amplitude');
    title('RESAMPLED TEMPLATE signals');
    grid on;
end
