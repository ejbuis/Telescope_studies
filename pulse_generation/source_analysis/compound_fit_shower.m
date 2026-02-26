% script that fits a combination of a cylindrical and a spherical source 
% amplitude decay profile to the 1e9 GeV neutrino shower decay profile
addpath('Datasets');
load("shower1e9.mat");
rpos = shower1e9.rpos;
p = shower1e9.p;
p = p/max(p);

splitind1 = 8;
splitind2 = 13;
f1 = fit(rpos(1:splitind1)', p(1:splitind1),'a/sqrt(x)', 'StartPoint', 1);
f2 = fit(rpos(splitind2:end)', p(splitind2:end),'a/x', 'StartPoint', 1);

r_fit = logspace(log10(min(rpos)), log10(max(rpos)), 100);
y_fit1 = f1(r_fit);
y_fit2 = f2(r_fit);

figure
hold on
plot(rpos, p, 'bx', 'DisplayName', '1e9 neutrino source')
plot(r_fit, y_fit1, 'b--', ...
    'DisplayName', '$\mathrm{Cylinder\ fit\ } \propto r^{-1/2}$')
plot(r_fit, y_fit2, 'r--', ...
    'DisplayName', '$\mathrm{Sphere\ fit\ } \propto r^{-1}$')

ylim([0 1])
xlabel('$r$ [m]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Peak pressure [a.u.]', 'Interpreter', 'latex', 'FontSize', 14)
set(gca, 'XScale', 'log', 'YScale', 'log', 'FontSize', 14);
legend('Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 12)

grid on
hold off