function rztp_tablefunc(Eo, gridvar)
    %% This script generates a binary table of (r,z,t,p)-values
    % To be used in Jpp implementations
    % Inputs:
    % Eo        [GeV]   Energy of the incoming neutrino
    % gridvar   [-]     Parameter for troubleshooting (determines gridsize)
    % toa       [-]     1: include time of arrival, 0: center all pulses around 0
    
    addpath('..');  
    addpath('../auxiliary');
    addpath('../source_analysis/Corsika_files');
    if nargin < 1; Eo = 1e9; end
    if nargin < 2; gridvar = 9; end
    
    c           = 1500;     % [m/s] Speed of sound
    fs          = 1e6;      % [Hz]  Sampling frequency
    wind        = 4096*2;   % 1/2   Size of the time window [-]
    savemat     = 1;        % [-]   Whether to save a mat file as well
    
    %% Spatial coordinates
    [rpos, zpos, t_axis, fn_spec] = generate_coordinates(gridvar, fs, wind);
    
    % Find required number of MC points from interpolation
    nmcs = nmc_from_r(rpos);
    print_array_info(rpos,zpos,t_axis,Eo);
    
    %% Initialize grid for storing data
    [robs, zobs] = meshgrid(zpos, rpos);
    [R, Z, T] = ndgrid(rpos, zpos, t_axis);
    ampgrid = zeros([length(rpos), length(zpos), length(t_axis)]);
    
    %% First generate the energy deposition:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fprintf("==== Generating From CORSIKA file ==== \n")
    % % => Binning specifics
    % cutdistance = 20; % [m]
    % [tsmc,Nz] = read_corsika_file('dEdzdr_51100_451.dat', cutdistance);
    % tsmc = tsmc / sum(sum(tsmc));
    % 
    % rsc = [.5:9.5 15:10:105];   % Radial bin centres (cm)
    % z_int = 20;                  % Interval of z (cm)
    % zsc = 10:20:(Nz*z_int);        % Longitudinal Bin Centres (cm)
    % fprintf('Total shower length %.1f\n', Nz*z_int/100);
    % fn_spec = strcat(fn_spec, '_reftail');
    % atten = 1;                % Learned's attenuation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    rsc = [.5:9.5 15:10:105]; % Radial bin centres (cm)
    zsc = 10:20:2000;         % Longitudinal Bin Centres (cm)
    atten = 7;                % ACoRNE attenuation - see atten_fna for the other options
    
    % Generate energy profile
    disp("NO MORE INTERPOLATING OF THE Eo GRID")
    tsmc = ShowerParm(rsc, zsc, Eo, 'CORSIKA');
    
    % NOTE: The following only works for the standard binning as included in
    % the Corsika parameterisation in Showerparm
    % tsmc = interpolate_Sm(Eo);
    
    % as the 10-100cm bins are 10x wider need to scale by a factor of 10 
    tsmc = tsmc * diag(kron([1 10], ones(1,10)));
    
    %% loop over the coordinates
    tic
    for j=1:length(zpos) 
    
        for i=1:length(rpos)
            Do = [rpos(i) 0 zpos(j)]; % Position of observer
    
            % Determine the monte carlo points, nmc = f(r)
            pointsc = MCGEn(tsmc, [0 zsc + 10], [0 rsc + [0.5 * ones(1, 10) 5 * ones(1, 10)]], nmcs(i));
    
            % Convert to cartesian
            [x_r, y_r, z_r] = pol2cart(rand(nmcs(i), 1) * 2 * pi, pointsc(:,2), pointsc(:,1));
    
            % Convert fom cm to m 
            points = [x_r y_r z_r] * 1e-2;
    
            % Compute the pressure
            [~,p,~,~] = kernelfrmod(points, Do, log10(Eo), atten, 10, fs, wind); % Perhaps think about generating just 1 instead of 10 rotations
            ampgrid(i,j,:) = p;
        end
        
        if rem(j, int8(length(zpos)/10)) == 0
            disp(strcat('Simulation at', 32,round(num2str(j/length(zpos)*100)), '%'));
        end
    
    end
    toc
    
    %% Plot
    amps_max = max(ampgrid, [], 3);
    figure()
    title(strcat('E = ', num2str(Eo/1e9), ' * 1e9 GeV'))
    surf(robs, zobs, amps_max)
    zscale('log')
    % xscale('log')
    xlabel('z [m]')
    ylabel('r [m]')
    zlabel('p_{max} [Pa]')
    
    %% Plot single trace
    r_ind = 2; z_ind = 2;
    t_a = sqrt((zpos(z_ind)-6)^2 + rpos(r_ind))^2/c;
    figure;
    hold on
    title(strcat('r = ', num2str(rpos(r_ind)), 'm, z = ', num2str(zpos(z_ind)), ...
        'm, toa = ', num2str(t_a), 's'))
    plot(squeeze(T(r_ind,z_ind,:)), squeeze(ampgrid(r_ind,z_ind,:))*1e6);
    xlabel('Time (s)');
    ylabel('Pressure (mPa)');
    % xlim([-100 100])
    % ylim([-12 12])
    grid on
    hold off
    
    
    %% Save data
    grid4D = cat(4,R,Z,T,ampgrid);
    tab = reshape(grid4D,[length(rpos)*length(zpos)*length(t_axis), 4]);
    tab_fix = reshape(tab.', [], 1);
    
    % File name generating
    fn = generate_filename(Eo, fn_spec);
    disp("Writing data to file")
    % Write to file
    fileID = fopen(fn, 'wb');
    fwrite(fileID, tab_fix, 'double');  % Save as double
    fclose(fileID);
    
    % Save in matlab format if desired
    if savemat
        fnmap = fn(1:end-4);
        fnmap = strcat(fnmap, '_ACORNEatten.mat');
        mapdata = struct();
        mapdata.r = rpos;
        mapdata.z = zpos;
        mapdata.t = T;
        mapdata.p = ampgrid;
        save(fnmap, 'mapdata')
    end
end