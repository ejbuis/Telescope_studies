function tsmc=disk_and_cylinder_source(r,z,w, z_maxs, sigma)
% Function that returns an energy distribution that consists of a cilinder
% along the total length of the z-binning + a gaussian "disk" at a provided peak
% value of z. The profile in r is Gaussian everywhere
% 
% INPUTS
% r         [cm]    bin centers in r
% z         [cm]    bin centers in z
% w         [a.u.]  relative weight cylinder & disk
% z_maxs    [cm]    [] location of maximum energy depositions
% sigma     [cm]    standard deviation of the location of maximum energy deposition
%
% OUTPUTS
% tsmc      [-]     normalised energy distribution over bins

% NOTE THIS ONLY WORKS WELL IF THE BINNING IN Z IS CONSTANT

% Generate & normalise cylinder, gaussian profile around r = 0
tsmc_cyl = ones(size(z(:)))*exp(-(4*(r/100)/sqrt(2)).^2); % QUESTION: r/100? 
tsmc_cyl = tsmc_cyl/sum(tsmc_cyl(:));
tsmc_disk = zeros(size(tsmc_cyl));

for i = 1:length(z_maxs)
    z_max = z_maxs(i);
    
    tsmci = w * exp(-(z - z_max).^2 / (2 * sigma^2))+1;
    tsmc_disk = tsmc_disk +  tsmci';

end

% Superimpose the energy maxima
tsmc = tsmc_cyl .* tsmc_disk;
tsmc = tsmc/sum(tsmc(:));