% Script to obtain a linear paramterization of the energy deposition maxima 
% for the shower parmaeterizations provided by ACORNE
addpath('..')
Eo              = logspace(8,12,10);
rscdiff         = [1,1,1,1,1,1,1,1,1,1,10,10,10,10,10,10,10,10,10,10];
rsc             = [.5:9.5 15:10:105];   % Radial bin centres (cm)
zsc             = 10:20:2000;           % Longitudinal Bin Centres (cm)
[rplot, zplot]  = meshgrid(rsc/100,zsc/100);

plotbool = 0;
zmaxs = zeros(size(Eo));

%%
for i = 1:length(Eo)
    
    tsmc = ShowerParm(rsc, zsc, Eo(i), 'CORSIKA');

    tsmcmod = zeros(size(tsmc));

    for j = 1:size(tsmc,1)
        tsmcmod(j,:) = tsmc(j,:) ./ rscdiff;
    end

    if plotbool
        figure;
        surf(rplot, zplot, tsmcmod);
    end

    [M,I] = max(tsmcmod(:,1));
    zmax = zsc(I);
    zmaxs(i) = zmax;
    
    % disp(strcat("Max energy deposition at z = ", num2str(zmax/100), " m"));

end

% Linear fit to the maximum energy deposition
a = polyfit(log10(Eo), zmaxs/100,1);
fprintf("fit parameters: a0 = %2.5f, a1 = %2.5f \n",a(2), a(1))

%%
figure
hold on
plot(Eo, zmaxs/100, 'x')
plot(Eo, polyval(a, log10(Eo)))
grid on
xscale('log')
xlabel('E0 [eV]')
ylabel('z_{max} [m]')