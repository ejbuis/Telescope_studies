function [rpos,zpos,t_axis, fn_spec] = generate_coordinates(gridvar, fs, wind)
% Function to generate a grid of points for pulse tabulating
%
% grid_var              What type of grid to use                [# points]
% 
%           1           Total grid (suitable for all energies)  [875]
%           2           Small grid                              [12]
%           3           Test grid                               [4]
%           4           Close to shower axis                    [63]
%           5           Rough sampling on a large grid          [49]
%
%           9           Optimized for 1e9 GeV                   [550]
%           10          Optimized for 1e10 GeV                  [650]                  
%           11          Optimized for 1e11 GeV                  [700]                  
%           12          Optimized for 1e12 GeV                  [700] 
%          else         Generic sparse grid                     [24]

%
% fs                    Sampling frequency               [Hz]
% wind                  1/2 size of the time window      [-]


% Generic grids
if gridvar == 1
    z1 = linspace(-100,-5,10);
    z2 = linspace(0,20,15);
    z3 = linspace(25,150,10);
    zpos = [z1,z2,z3];
    rpos = logspace(2,4,25); % now running from 100, to 10 km
    fn_spec = '';
elseif gridvar == 2
    zpos = [-10, 6, 20, 30];
    rpos = [100,1000,2000];
    fn_spec = '_small';
elseif gridvar == 3
    zpos = [0, 6];
    rpos = [1e2,1e3];
    fn_spec = '_xsmall';
elseif gridvar == 4
    zpos = [-100,-50,-25,-10,0,6,10,15,25,50,100];
    rpos = [0,1,5,10,25,40,50,75,100,125,150,175,200];
    fn_spec = '_closetoshower';
elseif gridvar == 5
    zpos = [-10000,-9000,-8000,-7000,-6000-100,0,100,5000,6000,7000,8000,9000,10000];
    rpos = [10,100,1000,10000,50000,100000,150000,200000,220000,250000,300000];
    fn_spec = '_roughlargegrid';

% Energy dependent grids
elseif gridvar == 9
    zpos = exponential_around_z6();
    rpos = [0,1,5,10,15,20,25,35,50,75,100,150,200,300,400,500,750,1e3,2e3];
    fn_spec = '_FULL';
    
elseif gridvar == 10
    zpos = exponential_around_z6();
    rpos = [0,1,5,10,15,20,25,35,50,75,100,150,200,300,400,500,600,750,850,1e3,1250,1500,2e3,3e3,5e3];
    fn_spec = '_FULL';

elseif gridvar == 11
    zpos = exponential_around_z6();
    rpos = [0,1,5,10,15,20,25,35,50,75,100,150,200,300,400,500,600,750,850,1e3,1250,1500,1750,2e3,2500,3e3,5e3,10e3];
    fn_spec = '_FULL';

elseif gridvar == 1200
    zpos = exponential_around_z6(50, -500, 500);
    rpos = [0,1,5,10,15,20,25,35,50,75,100,150,200,300,400,500,600,750,850,1e3,1250,1500,1750,2e3,2500,3e3,4e3,5e3,7500,10e3];
    fn_spec = '_FULL';

elseif gridvar == 12
    zpos = exponential_around_z6(50, -1e3, 1e3);
    rpos = [6000, 8000];
    fn_spec = '_1e12cleanup';

% Course end grid 
else
    zpos = [-6,0,6,12,20,50];
    rpos = [5e1, 1e2, 5e2, 1e3];
    fn_spec = '_rsparse';
end

% Time coordinates
t_axis = (-wind:wind-1)/fs; % time axis (centered at 0)
