function [x,y,z] = cylinder_uniform(nmc, z_l, sigma)
% Function that returns an energy distribution that consists of a cylinder
% running from -z_l/2 to z_l/2 centered in x and y with radius sigma.
% 
% INPUTS
% nmc       [-]     # of MC points to generate
% z_l       [cm]    length of cylinder in cm
% sigma     [cm]    radius of the cylinder
%
% OUTPUTS
% [x,y,z]   [cm]     coordinates of the points.

theta = 2*pi*rand(nmc,1);
radii = sigma*(rand(nmc,1).^(1/2));
x = radii .* cos(theta);
y = radii .* sin(theta);
z = z_l * (rand(nmc, 1) - 0.5);  

end
