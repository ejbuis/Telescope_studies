function [fn] = generate_filename(Eo, fn_spec)
% Function to generate a filename for pulse tabulating
%
% Eo                    Neutrino energy         [GeV]
% fn_spec               Gridsize specifier      [string]

fn = strcat('./rztpmaps_ACORNEatten/rztptab_', regexprep(sprintf('%.0e', Eo), 'e\+?0*', 'e'), 'GeV', fn_spec);

if exist(fn, 'file') == 2
    fnext = '_1';

    if exist(strcat(fn,fnext), 'file') == 2
        fnext = '_2';

        if exist(strcat(fn,fnext), 'file') == 2
            disp("Get your shit together");
            return
        end
    end

    fn = strcat(fn,fnext);
end

fn = strcat(fn, '.bin');
