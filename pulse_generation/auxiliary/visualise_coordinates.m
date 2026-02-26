function visualise_coordinates(rpos, zpos, fn_spec)
% VISUALISE_COORDINATES Scatterplot of all (r,z) coordinate combinations.
%
%   visualise_coordinates(rpos, zpos, fn_spec) creates a scatter plot of
%   all combinations of rpos and zpos (Cartesian product). The title
%   includes fn_spec to explain the data context.

    % Create full grid of coordinates using meshgrid
    [R, Z] = meshgrid(rpos, zpos);

    % Reshape into column vectors for plotting
    r_vals = R(:);
    z_vals = Z(:);

    % Create scatter plot
    figure;
    scatter(r_vals, z_vals, 'filled');
    grid on;
    xlabel('Radial position r');
    ylabel('Axial position z');
    xscale('log')

    % Title with explanation
    title(['Scatterplot of all (r,z) combinations - ' fn_spec]);
end
