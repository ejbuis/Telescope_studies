function visualise_rzprofile(rztpmap, p_thres)
% VISUALISE_RZPROFILE Visualizes the maximum pressure profile 
% in r-z coordinates and adds a transparent plane at p = 1e-3
%
% Inputs:
%   rztpmap - string: path to .mat file containing 'mapdata' struct
%   p_thres - scalar: threshold pressure for plotting horizontal plane

if nargin < 2
    plotpthres = 0;
else
    plotpthres = 1;
end

data = load(rztpmap);
if ~isfield(data, 'mapdata')
    error('The provided .mat file does not contain a "mapdata" struct.');
end

mapdata = data.mapdata;
title_red = split(rztpmap, '/');
nmred = title_red{end};
E0 = split(nmred, [".", "_"]);
E0 = E0{2};

% Extract fields
r = mapdata.r;  % 1xN
z = mapdata.z;  % 1xM
p = mapdata.p;  % NxMxT
[robs, zobs] = meshgrid(z, r);
amps_max = max(p, [], 3);

% Add transparent horizontal plane at p = 1e-3
if plotpthres
    z_plane = p_thres;
    [Z_plane_X, Z_plane_Y] = meshgrid(z, r);
    Z_plane_Z = ones(size(Z_plane_X)) * z_plane;
end

%% Plot
figure()
hold on
surf(robs, zobs, amps_max);

zscale('log');

if plotpthres
    h_plane = surf(Z_plane_X, Z_plane_Y, Z_plane_Z);
    set(h_plane, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [0.8 0.2 0.2]);
    legend('p_{max}', strcat('p = ', num2str(p_thres*1e3), ' mPa plane'), 'FontSize', 16);
else
    legend('p_{max}', 'FontSize', 20);
end

ylim([2 r(end)]); % r(end)
title(strcat('E = ', E0));
% fontsize(20, 'points');
xlabel('z [m]');
ylabel('r [m]');
xlim([-50  50]);
zlabel('p_{max} [Pa]');      % <-- add this line
view(3);         
caxis([1e-4 1e1]);
% cb = colorbar;
grid on;
set(gca, 'FontSize', 30, 'FontName', 'Times New Roman');
set(groot, 'defaultLineLineWidth', 2);
set(gca,'ColorScale','log')

hold off
end
