function tsmc = interpolate_Sm(E1)
    % Check E1 range
    if E1 < 1e5 || E1 > 1e12
        error('E1 must be between 1e5 and 1e12.');
    end
    
    addpath('..');
    % Load energy profile
    s = load('Sm');  % assumes s.Sm exists in Sm.mat

    % Convert E1 to index space
    E1ind = 2 * log10(E1) - 9;

    % Original index range
    rngE = 1:15;

    % Flatten s.Sm to 2D for interpolation
    [N, M, ~] = size(s.Sm);
    Sm_flat = reshape(s.Sm, [], 15);  % (N*M) x 15

    % Interpolate along the 3rd dimension (energy)
    Sm_interp_flat = interp1(rngE, Sm_flat.', E1ind, 'linear').';  % (N*M) x 1

    % Reshape back to NxM grid and transpose
    Sm_interp = reshape(Sm_interp_flat, N, M);
    tsmc = Sm_interp';  % return the transposed result
end
