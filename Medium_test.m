%% Testing K-wave to learn about Medium properties
% date: 21/10/2019
% Author: Jon Skerlj

clear all;
% Defining k-grid
Nx = 128;
dx = 0.1e-3; % grid point spacing in x dir [m/s]
Ny = 128;
dy = 0.1e-3; % grid point spacing in y dir [m/s]

% Making 2D grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% Defining a medium
c0 = 1000; % sound speed 1540 m/s
H20_den = 1000; % water density - kg/m^3

medium.sound_speed = c0 * ones(Nx,Ny); % making uniform sound speed for whole map
medium.sound_speed(:, 1:Ny/2) = 2000; % changing one part of the medium to other value

medium.density = 1000 * ones(Nx,Ny);
medium.density(:, 1:Ny/2) = 500;

medium.alpha_coeff = 0.75;  % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;

% create the time array
kgrid.makeTime(medium.sound_speed);


% create initial pressure distribution using makeDisc
% disc_magnitude = 6; % [Pa]
% disc_x_pos = 64;    % [grid points]
% disc_y_pos = 1;    % [grid points]
% disc_radius = 8;    % [grid points]
% disc_1 = disc_magnitude * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);
% 
% source.p0 = disc_1;

% % define a single source point
% source.p_mask = zeros(Nx, Ny);
% source.p_mask(end - Nx/4, Ny/2) = 1;
% 
% % define a time varying sinusoidal source
% source_freq = 0.25e6;   % [Hz]
% source_mag = 2;         % [Pa]
% source.p = source_mag * sin(2 * pi * source_freq * kgrid.t_array);

% define a single source point
source.p_mask = zeros(Nx, Ny);
source.p_mask(Nx/2, Ny/4) = 1;

% define a time varying sinusoidal source
source_freq = 0.25e6;   % [Hz]
source_mag = 2;         % [Pa]
source.p = source_mag * sin(2 * pi * source_freq * kgrid.t_array);

% filter the source to remove high frequencies not supported by the grid
source.p = filterTimeSeries(kgrid, medium, source.p);

% define a centered circular sensor
sensor_radius = 4e-3;   % [m]
num_sensor_points = 5;
sensor.mask = makeCartCircle(sensor_radius, num_sensor_points);

% define a 2D binary sensor mask in the shape of a line
% x_offset = 25;
% width = 50;
% sensor.mask = zeros();
% sensor.mask(x_offset, 60:70);

sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
    'PlotLayout', true, 'PlotPML', false);

% plot the simulated sensor data
figure;
imagesc(sensor_data, [-1, 1]);
colormap(getColorMap);
ylabel('Sensor Position');
xlabel('Time Step');
colorbar;