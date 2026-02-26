function [nmc] = nmc_from_r(r)
% Function that determines the number of MC points required at given
% distance(s) r. -> Thanks ChatGPT
%
% r         [m] observer location(s), can be scalar or array
%
% Reference points:
%   - At 100 m:      ~1e7 MC points -> IN FACT 1e8 is needed (unfeasible)
%   - At 1000 m:     ~1e6 MC points
%   - At 10000 m:    ~1e4 MC points

% Reference distances and corresponding MC points
r1 = 100;     nmc1 = 1e7;
r2 = 1000;    nmc2 = 1e6;
r3 = 10000;   nmc3 = 1e4;

% Preallocate output array
nmc = zeros(size(r));

% Region 1: r < r1
idx1 = r < r1;
nmc(idx1) = nmc1;

% Region 2: r >= r3
idx3 = r > r3;
nmc(idx3) = nmc3;

% Region 3: r1 <= r <= r2
idx12 = (r >= r1) & (r <= r2);
log_nmc12 = log10(nmc1) + (log10(nmc2) - log10(nmc1)) .* ...
            (log10(r(idx12)) - log10(r1)) ./ (log10(r2) - log10(r1));
nmc(idx12) = round(10.^log_nmc12, -3);

% Region 4: r2 < r <= r3
idx23 = (r > r2) & (r <= r3);
log_nmc23 = log10(nmc2) + (log10(nmc3) - log10(nmc2)) .* ...
            (log10(r(idx23)) - log10(r2)) ./ (log10(r3) - log10(r2));
nmc(idx23) = round(10.^log_nmc23, -3);

end
