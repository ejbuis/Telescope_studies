function visualise_isobar(rztpmap, p_iso)
    % VISUALISE_ISOBAR Visualizes the isobar surface at a given pressure
    % in 3D Cartesian coordinates with normalized axes for better visual balance.
    %
    % Inputs:
    %   rztpmap - string: path to .mat file containing 'mapdata' struct
    %   p_iso   - scalar: isobar pressure value to visualize

    % Load data
    data = load(rztpmap);
    if ~isfield(data, 'mapdata')
        error('The provided .mat file does not contain a "mapdata" struct.');
    end
    mapdata = data.mapdata;
    E0 = split(rztpmap, [".", "_"]);
    E0 = E0(end-1);


    % Extract fields
    r = mapdata.r;  % 1xN
    z = mapdata.z;  % 1xM
    p = mapdata.p;  % NxMxT

    % Get maximum pressure at each (r,z)
    p_max = max(p, [], 3);  % NxM

    [R, Z] = meshgrid(r, z);
    R = R'; Z = Z';

    % Define angular resolution
    nTheta = 100;
    theta = linspace(0, 2*pi, nTheta);

    % Convert normalized cylindrical to 3D Cartesian
    [ThetaGrid, RGrid] = meshgrid(theta, r);
    X = RGrid .* cos(ThetaGrid);  % (N x nTheta)
    Y = RGrid .* sin(ThetaGrid);  % (N x nTheta)

    % Build 3D grid
    numZ = length(z);
    X3D = zeros(length(r), nTheta, numZ);
    Y3D = zeros(length(r), nTheta, numZ);
    Z3D = zeros(length(r), nTheta, numZ);
    P3D = zeros(length(r), nTheta, numZ);

    for k = 1:numZ
        X3D(:, :, k) = X;
        Y3D(:, :, k) = Y;
        Z3D(:, :, k) = z(k);
        P3D(:, :, k) = repmat(p_max(:, k), 1, nTheta);
    end

    % Plot isosurface
    figure;
    pSurf = patch(isosurface(X3D, Y3D, Z3D, P3D, p_iso));
    set(pSurf, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.6);
    camlight; lighting gouraud;

    % Axis labels with physical units (but plotted in normalized space)
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    title(strcat('Isobar P = ', num2str(p_iso*1e3), 'mPa, E0 = ',E0, ' GeV'));

    fontsize(20, 'points');
    % axis equal;
    view(3); 
    grid on;
end
