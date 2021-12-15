clear

% Integrating mSOUND with k-Wave for transducers of curved shape
%% ======================================================================== 
%%                        PARAMETERS 
%% ========================================================================

fc = 1.0e6;      % fundamental frequency [Hz]
c0 = 1500;       % speed of Sound [m/s]
p0 = 1;   
rho0 = 1000;     % density of medium [kg/m^3]
lambda = c0/fc;  % wavelength [m]

dy = lambda/8;   % step size in the x direction [m]
dx = lambda/8;   % step size in the x direction [m]

% assign the grid size and create the computational grid
Nx = 240;   
% Ny_p = 80;
Ny_p = 30;    % actural propagation distance
Ny = Ny_p + 20;

kgrid = kWaveGrid(Nx, dx, Ny, dy);
%% create the time array
dt = 1/fc/32;   % time step in kwave [s]
% cfl = c0*dt/dx;

kgrid.t_array = 0:dt:4/fc*10;  % time array [s]

%% assign the properties of the propagation medium
medium.sound_speed = 1500;
medium.BonA        = 2.0*(0.0-1.0);
medium.density     = 1000;
medium.alpha_coeff = 0.0;
medium.alpha_power = 2.0;

%% define the source and the sensor
source.p_mask = zeros(Nx, Ny);
sensor.mask = zeros(Nx, Ny);
offcenter = floor(Ny/2 - 10);

TR_radius = 8*lambda;         % transducer radius [m]  
TR_focus  = 10*lambda;        % transducer focal length [m]
focus  = round(TR_focus/dy);  % [grid points]

cx = round(Nx/2)+1;
cy = 10 + focus;

arc_angle1 = asin(TR_radius/TR_focus);
circle1 = makeCircle(Nx, Ny, cx, cy, focus, arc_angle1, false);
arc_angle2 = pi*2;
circle2 = makeCircle(Nx, Ny, cx, cy, focus, arc_angle2, false);
arc_angle3 = pi*2 - arc_angle1;
circle3 = makeCircle(Nx, Ny, cx, cy, focus, arc_angle3, false);
circle = circle2 - circle3 + circle1;

source.p_mask = circle;
sensor.mask(:, Ny-10) = 1; 

figure
imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, circle)
axis image
%% excitation signal

ts = [-4/fc:dt:4/fc];
excit_p = sin(2*pi*fc*(ts)).*exp(-(ts).^2*fc^2/2);
source.p = excit_p*p0;

input_args  = {'DisplayMask', source.p_mask, 'DataCast', 'single'};
% run the simulation with kwave
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
%%
sensor_data = reshape(sensor_data,  Nx,length(kgrid.t_array));
figure
plot(kgrid.t_array*1e6, squeeze(sensor_data( round(Nx/2), :)), 'r', 'LineWidth', 2)

%% ======================================================================== 
%%                        Connect to MDM
%% ========================================================================
x_length = kgrid.dx*kgrid.Nx;
y_length = TR_focus - kgrid.dy*(Ny-20);
t_length = max(kgrid.t_array); 

dt = 1/fc/16;  % timestep in TMDM [s]

mgrid = set_grid(dt, t_length, dx, x_length, dy, y_length);

source_p = sensor_data(:, 1:2:end).';

medium.c    = 1500;
medium.rho  = 1000;
medium.beta = 0.0;
medium.ca   = 0;       % power law relation constant     
medium.cb   = 2.0;

medium.c0 = 1500;      % reference speed of Sound [m/s]

medium.NRL_gamma = 0.5;
medium.NRL_alpha = 0.05;

% define where the pressure will be recorded
sensor_mask = zeros(mgrid.num_x, mgrid.num_y+1);
sensor_mask(:, end) = 1;

% 2D forward simulation
p_time = Forward2D(mgrid, medium, source_p, sensor_mask, 0);

%% ======================================================================== 
%%                        RESULTS
%% ========================================================================
p_time = reshape(p_time, mgrid.num_t, mgrid.num_x);

figure
p_focus(:) = p_time(:, round(mgrid.num_x/2));
plot(mgrid.t*1e6, p_focus, 'b', 'LineWidth',2)
xlabel ('Time (\mus)')
ylabel ('Pressure (Pa)')





