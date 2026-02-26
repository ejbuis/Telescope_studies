% Script to resample a list rztp maps

allEs = {...
    '../rztpmaps_ACORNEatten/rztptab_1e9GeV_ACORNEatten.mat', ...
    '../rztpmaps_ACORNEatten/rztptab_5e9GeV_ACORNEatten.mat', ...
    };

%% Downsample

f0 = 1e6;     % original sampling rate (Hz)
f1 = 144e3;   % target sampling rate (Hz)

Nt_all = zeros(numel(allEs),1);

for i = 1:numel(allEs)
    fprintf('\n--- Processing %d / %d ---\n', i, numel(allEs));
    Nt_all(i) = resample_rztpmap(allEs{i}, f0, f1,false, false);
end

disp('All files processed successfully!');
disp(table(allEs', Nt_all, 'VariableNames', {'Filename','Nt'}));
