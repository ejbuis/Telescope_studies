function [b,a] = hydrophone(G_dB, G_res_dB, f_res, fw_res, fs_filt)
%HYDROPHONE: returns transfer function consisting of a constant gain and a
% resonance peak, as a first order approximation of the hydrophone
% characteristic.
%
% H(z) = B(z)         gr(1-R^2)                 1-z^-2
%        ----  = g0 + --------- * --------------------------------
%        A(z)             2        1 - 2Rcos(theta)z^-1 + R^2 z^-2
% 
%  (https://ccrma.stanford.edu/~jos/filters/Constant_Peak_Gain_Resonator.html)
% 
% INPUTS:
% G_dB      [dB]: constant gain
% G_res_dB  [dB]: resonance peak gain
% f_res     [Hz]: frequency of the resonance peak
% fw_res    [Hz]: FWHM of the resonance peak
% fs_filt   [Hz]: sampling frequency of the filter

if nargin < 5 fs_filt = 1e6; end

% Convert gains from dB to linear
g               = 10^(G_dB / 20);       % Baseline gain
g_res           = 10^(G_res_dB / 20);   % Total peak gain at resonance
g_r             = g_res - g;            % Peak amplitude *above* baseline

% Compute radius and angle
R       = exp(-pi*fw_res/fs_filt);
theta   = 2*pi*f_res/fs_filt;

% compute filter coefficients
b       = [1, 0, -1];
a       = zeros(1,3);
a(1)    = 1;
a(2)    = -2 * R * cos(theta);
a(3)    = R^2;

% scale peak gain to 1
b = b .* (1-R^2)/2;

% apply relative peak gain
b = b .* g_r;

% apply base gain
b = b + g*a;

end

