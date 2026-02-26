function print_array_info(rpos, zpos, t_axis, Eo)
% Prints formatted information about r, z, t arrays and E₀ in GeV.
fprintf('\n======================\n')
fprintf('Writing (rztp)-table for:\n');

% Get array lengths
nr = length(rpos);
nz = length(zpos);
nt = length(t_axis);

% Get min/max values
rmin = min(rpos); rmax = max(rpos);
zmin = min(zpos); zmax = max(zpos);

% Print sizes and ranges
fprintf('  r (%4d) : %.2f m to %.2f m\n', nr, rmin, rmax);
fprintf('  z (%4d) : %.2f m to %.2f m\n', nz, zmin, zmax);
fprintf('  t (%4d) : %.2e s to %.2e s\n', nt, min(t_axis), max(t_axis));

% Convert energy to GeV and round
Eo_GeV_rounded = round(Eo, 3 - floor(log10(abs(Eo))) - 1);
fprintf('  E0       : %.3g GeV\n', Eo_GeV_rounded);
fprintf('======================\n')
end
