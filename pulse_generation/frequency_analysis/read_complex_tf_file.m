fid = fopen('./Transfer_functions/Datamat_Dubbel_1FF_Atran_Sweep200Hz20kHz180degreev3.dat', 'r');

fs = [];
H_complex = [];

pattern = '\(\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s*([+-])\s*([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)j\s*\)\s+\(\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s*([+-])\s*([0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)j\s*\)';

while ~feof(fid)
    line = fgetl(fid);
    tokens = regexp(line, pattern, 'tokens');
    
    if ~isempty(tokens)
        t = tokens{1};  % cell array of strings
        
        % Convert freq real and imag
        freq_real = str2double(t{1});
        freq_sign = t{2};
        freq_imag = str2double(t{3});
        if freq_sign == '-'
            freq_imag = -freq_imag;
        end
        
        % Convert response real and imag
        resp_real = str2double(t{4});
        resp_sign = t{5};
        resp_imag = str2double(t{6});
        if resp_sign == '-'
            resp_imag = -resp_imag;
        end
        
        % Append
        fs(end+1,1) = complex(freq_real, freq_imag);
        H_complex(end+1,1) = complex(resp_real, resp_imag);
    end
end
fs = real(fs);
fclose(fid);
save("./Transfer_functions/reponse180deg","H_complex", "fs")