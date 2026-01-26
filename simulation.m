clearvars; close all; clc;

source_freq = 55e3;    % time varying source
source_mag = 1e4;
c_source = 1500;          % water
lambda = c_source / source_freq;

%% creating computational grid
Nx = 216;
Ny = 216;
points_per_wavelength = 10;
dx = lambda/points_per_wavelength; % must resolve shortest wavelength, less than lambda/10
dy = dx;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

%% creating homogenous medium (water)
medium.sound_speed = 1500 * ones(Nx, Ny);
medium.density = 1000 * ones(Nx, Ny);

%% create circular target
target_radius = 10;
target_x = round(Nx * .8);  % towards the bottom of the grid
target_y = round(Ny / 2);   % centered
target_mask = makeDisc(Nx, Ny, target_x, target_y, target_radius);

% impedance mismatch for reflections (air)
medium.sound_speed(target_mask == 1) = 343;
medium.density(target_mask == 1) = 1.2;

%% creating metal walls
margin = ceil(lambda/(4*dx));
wall_thickness = 5;

% bottom wall
bottom_wall = false(Nx, Ny);
bottom_wall(end-wall_thickness+1-margin:end-margin, :) = true;

% left wall
left_wall = false(Nx, Ny);
left_wall(:, margin+1:wall_thickness+margin) = true;

% right wall
right_wall = false(Nx, Ny);
right_wall(:, end-wall_thickness+1-margin:end-margin) = true;

% combine all walls
tank_walls = bottom_wall | left_wall | right_wall;

% apply steel properties to walls
medium.sound_speed(tank_walls) = 5900;
medium.density(tank_walls) = 7900;

% %% creating dampener
% 
% medium.alpha_coeff = zeros(Nx, Ny);  % initialize attenuation for dampener
% medium.alpha_power = 1.5;
% damp_thickness = 8;  % thicker than steel wall
% 
% % bottom dampener
% bottom_dampener = false(Nx, Ny);
% bottom_dampener(end-wall_thickness-damp_thickness-margin+1:end-wall_thickness-margin, :) = true;
% 
% % left dampener
% left_dampener = false(Nx, Ny);
% left_dampener(:, wall_thickness+margin+1:wall_thickness+damp_thickness+margin) = true;
% 
% % right dampener
% right_dampener = false(Nx, Ny);
% right_dampener(:, end-wall_thickness-damp_thickness-margin+1:end-wall_thickness-margin) = true;
% 
% % combine dampeners
% dampeners = bottom_dampener | left_dampener | right_dampener;
% 
% % dampener properties (rubber/foam-like)
% medium.sound_speed(dampeners) = 1800;      
% medium.density(dampeners) = 1200;          
% medium.alpha_coeff(dampeners) = 5.0;       % HIGH attenuation (dB/(MHz^y cm))

%% time array
c_max = max(medium.sound_speed(:));
t_end = 2 * Nx * dx / c_source;   % enough time for propogation 2 directions
kgrid.t_array = makeTime(kgrid, c_max, 0.3, t_end); % CFL = 0.3 

%% point source
buffer = 20 + ceil(lambda / (4*dx));

source.p_mask = zeros(Nx, Ny);
source_x = buffer;
source_y = round(Ny * 0.3);
source.p_mask(source_x, source_y) = 1;

source.p = source_mag * toneBurst(1/kgrid.dt, source_freq, 5);

%% point sensor
sensor.mask = zeros(Nx, Ny);
sensor_x = buffer; % at least lambda/4 from PML
sensor_y = round(Ny * .7);
sensor.mask(sensor_x, sensor_y) = 1;  
sensor.record = {'p'};

%% see tank visualization
figure;
imagesc(medium.density);
axis image;
colorbar;
title('Simulation Setup');
hold on;

% source
[src_y, src_x] = find(source.p_mask);   % row, col
plot(src_x, src_y, 'r*', 'MarkerSize', 10, 'LineWidth', 1.5);

% sensor
[sen_y, sen_x] = find(sensor.mask);
plot(sen_x, sen_y, 'bo', 'MarkerSize', 8, 'LineWidth', 1.5);
legend('Source', 'Sensor');
hold off;

%% running simulation
cm = getColorMap();
input_args = {'PMLInside', false, 'PMLSize', 20, 'ColorMap', cm, 'DataCast', 'single'};
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

%% plot A-scan of point sensor
figure;
plot(kgrid.t_array*1e6, sensor_data.p, 'k');  % time in microseconds
xlabel('Time [\mus]');
ylabel('Pressure [Pa]');
title('A-scan at Single Sensor');
grid on;