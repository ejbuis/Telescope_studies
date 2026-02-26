function zvals = exponential_around_z6(N, z_min, z_max)
    
    if nargin < 1; N = 50; end
    if nargin < 2; z_min = -200; end
    if nargin < 3; z_max = 200; end

    mu      = 6;              % Gaussian center
    lambda  = 5;                % exponential decay rate (bigger = stronger clustering)
    
    % Number of points left and right of the center (excluding center itself)
    nL = floor((N-1)/2);
    nR = N - nL - 1;
    
    % Create (nL+1) samples from mu -> x_min (first is mu, last is x_min)
    tL = linspace(0,1,nL+1);
    dL = (exp(lambda * tL) - 1) / (exp(lambda) - 1);   % 0..1
    xL_full = mu - dL * (mu - z_min);                  % length nL+1: [mu ... x_min]
    
    % Take elements 2:end and reverse to get ascending order from x_min up to (left of mu)
    xL = fliplr(xL_full(2:end));                       % length nL
    
    % Create (nR+1) samples from mu -> x_max (first is mu, last is x_max)
    tR = linspace(0,1,nR+1);
    dR = (exp(lambda * tR) - 1) / (exp(lambda) - 1);   % 0..1
    xR_full = mu + dR * (z_max - mu);                  % length nR+1: [mu ... x_max]
    
    % Take elements 2:end (already ascending from just-right-of-mu -> x_max)
    xR = xR_full(2:end);                               % length nR
    
    % Combine: left ascending, center, right ascending
    zvals = [xL, mu, xR];                                 % length nL + 1 + nR == N

end
