function [t_return, p,pw,Exyz]=kernelfrmod(points,Do,energy,atten,nr,fs,wind, m,c)
% UltraFast Acoustic Integral function -> modified 20250424 by Huib Baetsen 
% + readability, + return uncorrected time axis, + time window as variable
%
% INPUTS
% points:   [x,y,z] where z is oriented along the axis of the shower,
%           units (m) size  n x 3 where n is typically about 10^6. Note the value of n is very
%           distance dependent. At 100m c 10^8 points are needed to give an ultra-clean pulse.
%           At 10km c 10^4 points are sufficient.
%
% Do:       [x0,y0,z0] the position of the observer 
%
% energy    total energy of the shower 10^energy GeV. i.e. if energy = 11
%           then the total shower energy is 10^20 eV.
%
% atten     is a flag and is passed to attenfna.
%
% nr        rotational symmetry is exploited by rotating the shower axially. a
%           value of 100 is typical. (default 1)
% 
% fs        sampling frequency (default 1MHz)
%
% wind      1/2 size of the time window (default 1024)
% 
% m         mean distance from observer to shower. Calculated if not supplied
% 
% OUTPUTS
% p         the pulse pressure (sampling rate 1MHz default) Note: zero time is at the mean shower
%           transit time ignoring complex attenuation
%
% pw        the FFT of the pulse (sampling rate 1MHz default)
%
% Exyz      a scaled version of the Velocity Potential
%
% SD Last mod 7/7/08

% ------------------------------------------------------------------------- 
%% DEFINE CONSTANTS
% Velocity of sound in water(c),thermal Expansivity (beta), Specific heat capacity (Cp)
beta = 2.0e-4;      
Cp = 3.8e3;
if nargin <5 nr = 1;    end
if nargin <6 fs = 1e6;  end 
if nargin <7 wind = 1024; end
if nargin <9 c = 1500;  end

t_axis = (-wind:wind-1)'/fs;
f_axis = [0:(wind-1) -wind:-1]'/(2*wind) * fs; %Time and frequency axes -> was /2028 ?!

% differential in Frequency Domain. The d/dt of exp(iwt) is iw exp(iwt)
% A Blackman window is used to smooth the integral and is optional
diff_filt = 1i * 2 * pi * f_axis .* fftshift(blackman(2*wind));

% Create Velocity potential array
Exyz = zeros(2*wind,1);

% -------------------------------------------------------------------------
% Determine number of points in the integral
nmc = size(points, 1);
% If rotational symmetry is being used calculates the angles
th = linspace(0,2*pi,1+nr);
th(end) = [];
% Convert the observer position to a matrix of the same size as the number
% of points
Do = Do(ones(1, nmc), :);

% loop through the angles
for ang = th
    % Rotate using matrix rotation 
    pointsr = points * [cos(ang)  sin(ang) 0; ...
                        -sin(ang) cos(ang) 0; ...
                        0         0        1];
    % determine the distance to each point
    d = sqrt(sum((pointsr - Do) .^ 2, 2));

    % Determine the mean distance if not provided
    if nargin<8
        m=sum(d)/nmc; 
    end

    % Deterime the unscaled velocity potential. (Need to divide by distance)
    his = histc(d/c, t_axis + m/c);
    Exyz = Exyz + his(:);
end
% figure
% hold on
% title("distance to MC points - mean distance")
% xlabel("d-m [m]")
% histogram(d-m,200)
% hold off

% Normalise sum(Exyz) to one. 
% histogram(Exyz, length(t_axis));
Exyz = Exyz / length(th) / nmc;

% Scale by constants and divide by distance. 
Exyzn = Exyz * beta / Cp / 4 / pi * 10^(9+energy) * 1.6e-19 ./ (m + t_axis*c);

% Do integral in the frequency domain
pw = fft(Exyzn) .* diff_filt .* atten_fna(f_axis, m*1e-3, atten);

% Convert back to the time domain
p = real(ifft(pw) * fs);
% disp(strcat("m = ", num2str(m), " m"));
t_return = t_axis + m/c;
