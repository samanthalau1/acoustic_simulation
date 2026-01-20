%% creating computational grid
Nx = 256;
Ny = 256;
dx = 1e-4
dy = dx;

kgrid = kWaveGrid(Nx, dx, Ny, dy);

%% creating homogenous medium (water)
medium.sound_speed = 1500 * ones(Nx, Ny);
medium.density = 1000 * ones(Nx, Ny);

%% create circular target
target_radius = 10;
target_x = round(Nx * .8);  % towards the bottom of the grid
target_y = round(Ny / 2);   % centered

target_mask = makeCircle(Nx, Ny, target_x, target_y, target_radius);

% impedance mismatch for reflections
medium.sound_speed(target_mask == 1) = 5900;
medium.density(target_mask == 1) = 7900;

%% creating metal walls
wall_thickness = 5;

% bottom wall
bottom_wall = false(Nx, Ny);
bottom_wall(end-wall_thickness+1:end, :) = true;

% left wall
left_wall = false(Nx, Ny);
left_wall(:, 1:wall_thickness) = true;

% right wall
right_wall = false(Nx, Ny);
right_wall(:, end-wall_thickness+1:end) = true;

% combine all walls
tank_walls = bottom_wall | left_wall | right_wall;

% apply steel properties to walls
medium.sound_speed(tank_walls) = 5900;
medium.density(tank_walls) = 7900;

%% time array
c_max = max(medium.sound_speed(:));
t_end = 2 * Nx * dx / 1500;
kgrid.makeTime(c_max, 0.3, t_end);

%% point source
source.p_mask = zeros(Nx, Ny);

source_x = 2;        % below surface
source_y = round(Ny * .3);  % left of sensor
source.p_mask(source_x, source_y) = 1;

source_freq = 200e3;   % time varying source
source_mag  = 5e5;

source.p = source_mag * toneBurst(1/kgrid.dt, source_freq, 1);

%% point sensor
sensor.mask = zeros(Nx, Ny);
sensor_x = 2;                       % very top = surface
sensor_y = round(Ny * .7);
sensor.mask(sensor_x, sensor_y) = 1;

sensor.record = {'p'};

%% see tank
figure;
imagesc(medium.density');
colorbar;
title('Tank Geometry (density)');

%% running simulation
input_args = {'PMLInside', false, 'PMLSize', 40, 'DataCast', 'single'};
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

%% plot A-scan of point sensor
figure;
plot(kgrid.t_array*1e6, sensor_data.p, 'k');  % time in microseconds
xlabel('Time [\mus]');
ylabel('Pressure [Pa]');
title('A-scan at Single Sensor');
grid on;

