function Nt_new = resample_rztpmap(filename, f0, f1,doRescale, doPlot)
%RESAMPLE_RZTPMAP Resample an r-z-t pressure map and save as binary file
%
%   Nt_new = RESAMPLE_RZTPMAP(filename, f0, f1, doPlot)
%
%   INPUTS:
%       filename : string
%           Path to the input .mat file containing mapdata
%       f0 : double
%           Original sampling frequency [Hz] (e.g. 1e6)
%       f1 : double
%           Target sampling frequency [Hz] (e.g. 144e3)
%       doRescale : logical (optional)
%           If true, the dataset is rescaled such that the max amplitude
%           at each gridpoint matches that of the f0 map.
%       doPlot : logical (optional)
%           If true, plots one example (r,z) trace before/after resampling
%           Default = false
%
%   OUTPUTS:
%       Nt_new : integer
%           Number of points in the resulting resampled time trace
%
%   The function loads the r-z-t pressure map, resamples the time traces
%   from f0 to f1, rescales the amplitudes to conserve peak pressure,
%   and saves the result into a binary file
%   "<filename>_resampled_<f1>kHz.bin".
%

if nargin < 2; f0        = 1e6;      end
if nargin < 3; f1        = 144e3;    end 
if nargin < 4; doRescale = false;    end
if nargin < 5; doPlot    = false;    end

% --- Check if file exists ---
assert(isfile(filename), 'Input file "%s" does not exist.', filename);

% --- Load data ---
dat         = load(filename);
mapdata     = dat.mapdata;
pds         = mapdata.p;
ts          = mapdata.t;
rs          = mapdata.r;
zs          = mapdata.z;
maxps_orig  = max(pds, [], 3);

% --- Resample entire map ---
pdsamp      = resample(pds, f1, f0, 2, 'Dimension', 3);
maxps_res   = max(pdsamp, [], 3);
dpsDownSamp = abs(maxps_res ./ maxps_orig);

fprintf('Downsampling / %.0f Hz: mean %.4f, stddev %.4f \n', ...
    f0, mean(dpsDownSamp,"all"), std(dpsDownSamp,0,"all"));

% --- Rescale to conserve amplitude ---
if doRescale
    rescaled_map        = pdsamp ./ dpsDownSamp;
    maxps_rescaled      = max(rescaled_map, [], 3);
    dpsAfterSamp        = abs(maxps_rescaled ./ maxps_orig);
    fprintf("Rescaling performed, new mean amplitude ratio = %.3f \n", mean(mean(dpsAfterSamp)));
else
    rescaled_map        = pdsamp;
end

% --- Plot one trace if requested ---
if doPlot
    r_ind   = 3; % example index in r
    z_ind   = 30; % example index in z
    p_orig  = squeeze(pds(r_ind, z_ind, :));
    p_rescl = squeeze(rescaled_map(r_ind, z_ind, :));
    p_resmp = squeeze(pdsamp(r_ind, z_ind, :));
    t_orig  = squeeze(ts(r_ind, z_ind, :));
    t_resmp = linspace(min(t_orig), max(t_orig), length(p_resmp));

    figure;
    hold on
    plot(t_orig,  p_orig,  '-o', 'DisplayName', 'Original signal')
    plot(t_resmp, p_rescl, '-o', 'DisplayName', 'Resampled, scaled')
    plot(t_resmp, p_resmp, '-o', 'DisplayName', 'Resampled raw')
    title(sprintf("Traces for r = %d m, z = %d m", rs(r_ind), zs(z_ind)))
    grid on
    legend
    xlabel("t [s]")
    ylabel("P [Pa]")
    set(gca, 'FontSize', 12)
    hold off
end

% --- Save the rescaled map in binary format ---
Nt_new      = size(rescaled_map, 3);
indices     = linspace(0, Nt_new-1, Nt_new);

[R, Z, I]   = ndgrid(rs, zs, indices);
grid4D      = cat(4, R, Z, I, rescaled_map);

tab         = reshape(grid4D, [length(rs)*length(zs)*Nt_new, 4]);
tab_fix     = reshape(tab.', [], 1);

% Format f1 in kHz for filename
f1_kHz = f1 / 1e3;
if abs(f1_kHz - round(f1_kHz)) < eps
    f1_str = sprintf('%d', round(f1_kHz));  % integer kHz
else
    f1_str = sprintf('%.3g', f1_kHz);       % scientific kHz
end

% --- Build output path ---
[~, name, ~] = fileparts(filename);

% Replace any trailing '_totalgrid_...' or similar with '_resampled_<f1>kHz'
outname_base = regexprep(name, '_totalgrid_.*$', sprintf('_resampled_%skHz', f1_str));

if doRescale
    outdir = './resampled_maps/rescaled';
    outname = fullfile(outdir, [outname_base '_rescaled.bin']);
else
    outdir = './resampled_maps';
    outname = fullfile(outdir, [outname_base '.bin']);
end

if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% Write binary
fileID = fopen(outname, 'wb');
fwrite(fileID, tab_fix, 'double');  % Save as double precision
fclose(fileID);

fprintf('Saved resampled map to %s\n', outname);

end
