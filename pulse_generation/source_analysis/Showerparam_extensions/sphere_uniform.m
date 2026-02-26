function [x,y,z] = sphere_uniform(nmc, z_max, sigma)
% Function that returns an energy distribution that consists of a spherical Gaussian
% centered at (r=0, z=z_max) with standard deviation sigma.
% 
% INPUTS
% r         [cm]    bin centers in r
% z         [cm]    bin centers in z
% ~         Ignored (kept for compatibility with original call signature)
% z_max     [cm]    location of maximum energy deposition
% sigma     [cm]    standard deviation of the energy deposition
%
% OUTPUTS
% tsmc      [-]     normalised energy distribution over bins

rvals = 2*rand(nmc,1)-1;
elevation = asin(rvals);
azimuth = 2*pi*rand(nmc,1);
radii = sigma*(rand(nmc,1).^(1/3));

[x,y,z] = sph2cart(azimuth,elevation,radii);
z = z+z_max;
end
