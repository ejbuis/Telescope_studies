function [fit_result, y_fit, params] = fit_exp_retention(f, y)
%FIT_EXP_RETENTION Fits exponential amplitude retention curve:
%   y(f) = 1 - C * exp(-f / fc)
%
% Inputs:
%   f - frequency data (vector)
%   y - normalized amplitude retention (vector)
%
% Outputs:
%   fit_result - MATLAB fit object
%   y_fit - fitted values at input frequencies
%   params - structure with fields:
%       C  - fit coefficient C
%       fc - fit coefficient fc (characteristic frequency)

    % Ensure column vectors
    f = f(:);
    y = y(:);

    % Define fit type
    fit_model = fittype('1 - C*exp(-x/fc)', ...
        'independent', 'x', ...
        'coefficients', {'C', 'fc'});

    % Fit options with bounds
    opts = fitoptions(fit_model);
    opts.StartPoint = [0.5, 3e5];       % Start with fc ~ 300 kHz
    opts.Lower = [0, 1e4];              % fc at least 10 kHz
    opts.Upper = [1.5, Inf];            % No upper bound

    % Fit the model
    fit_result = fit(f, y, fit_model, opts);

    % Evaluate fitted values
    y_fit = fit_result(f);

    % Output parameters
    if nargout > 2
        params.C = fit_result.C;
        params.fc = fit_result.fc;
    end
end
