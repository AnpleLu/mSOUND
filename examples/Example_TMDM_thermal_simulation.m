clear
% Integrating mSOUND with k-Wave for thermal simulations
%% ======================================================================== 
%%                        PARAMETERS 
%% ========================================================================
fc = 1.0e6;	           % fundamental frequency[Hz]
p0 = 0.6e6;	           % initial pressure amplitude [Pa]
medium.c0 = 1500;      % speed of Sound [m/s]
medium.rho0 = 1000;    % density of medium [kg/m^3]

omegac = fc*2*pi;      % fundamental angular frequency
lambda = medium.c0/fc; % wavelength [m]

dx = lambda/10;        % step size in the x direction [m]
dy = lambda/10;        % step size in the y direction [m]
dt = 1/fc/16;          % timestep [s]

x_length = lambda*30+dx; % computational domain size in the x direction [m]
y_length = lambda*20;    % computational domain size in the y direction [m]
t_length = 4/fc*15.0;    % temporal domain size [s] 
%% ======================================================================== 
%%                        COMPUTATIONAL DOMAIN
%% ========================================================================

mgrid = set_grid(dt, t_length, dx, x_length, dy, y_length);  %% wavevector and x, y

%% ======================================================================== 
%%                        TRANSDUCER
%% ========================================================================

TR_focus  = 12*lambda;  % transducer focal length [m]
TR_radius = 6*lambda;   % Transducer radius   [m]

%% ======================================================================== 
%%                        EXCITATION
%% ========================================================================

% calculate the phase delay
delay = sqrt((mgrid.x).^2 + (TR_focus)^2)/medium.c0;  % [s]
delay = delay - min(delay);

delay = repmat(delay, mgrid.num_t,1);
ts = repmat(mgrid.t.', 1,mgrid.num_x);

% define the excitation pulse
source_p = p0*sin(2*pi*fc*(ts+delay));
source_p(:, abs(mgrid.x)>TR_radius) = 0;
%% ======================================================================== 
%%                        MEDIUM
%% ========================================================================
medium.c    = 1500;   % speed of sound [m/s]      
medium.rho  = 1000;   % density [kg/m^3]             
medium.beta = 0;      % nonlinear coefficient     
medium.ca   = 0.75;   % attenuation coefficient [dB/(MHz^y cm)]     
medium.cb   = 1.5;    % power law exponent    

% non-reflecting layer
medium.NRL_gamma = 0.5;
medium.NRL_alpha = 0.1;
%% ======================================================================== 
%%                        SIMULATION
%% ========================================================================

% define where the pressure will be recorded
sensor_mask = zeros(mgrid.num_x, mgrid.num_y+1);
sensor_mask(:, 2:end) = 1;

% 2D forward simulation
p_time = Forward2D(mgrid, medium, source_p, sensor_mask, 0, 'NRL');

p = extractAmpPhase(p_time.', 1/dt, fc);

% reshape the data, and calculate the volume rate of heat deposition
p = reshape(p, mgrid.num_x, mgrid.num_y);
p = abs(p);

%% ======================================================================== 
%%                        thermal simulation
%% ========================================================================

% convert the absorption coefficient to nepers/
alpha_np = db2neper(medium.ca, medium.cb) * (omegac).^medium.cb;
Q = alpha_np .* p.^2 ./ (medium.rho .* medium.c);

% set the background temperature and heating term
source.Q = Q;
source.T0 = 37;

% clear medium
% define medium properties related to diffusion
medium.density              = 1000;     % [kg/m^3]
medium.thermal_conductivity = 0.5;      % [W/(m.K)]
medium.specific_heat        = 3600;     % [J/(kg.K)]

% create the computational grid
kgrid = kWaveGrid(mgrid.num_x, mgrid.dx, mgrid.num_y, mgrid.dy);

% create kWaveDiffusion object
kdiff = kWaveDiffusion(kgrid, medium, source, []);

% set source on time and off time
on_time  = 10;  % [s]
off_time = 20;  % [s]
% set time step size
dt = 0.1;

% take time steps
kdiff.takeTimeStep(round(on_time / dt), dt);

% store the current temperature field
T1 = kdiff.T;

% turn off heat source and take time steps
kdiff.Q = 0;
kdiff.takeTimeStep(round(off_time / dt), dt);

% store the current temperature field
T2 = kdiff.T;

%% ======================================================================== 
%%                        VISUALISATION
%% ========================================================================
% plot the temperature after heating
figure;

% plot the acoustic pressure
subplot(2, 3, 1);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, p * 1e-6);
h = colorbar;
xlabel(h, '[MPa]');
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Acoustic Pressure Amplitude');

% plot the volume rate of heat deposition
subplot(2, 3, 2);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, Q * 1e-7);
h = colorbar;
xlabel(h, '[kW/cm^3]');
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Volume Rate Of Heat Deposition');

% plot the temperature after heating
subplot(2, 3, 3);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, T1);
h = colorbar;
xlabel(h, '[degC]');
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Temperature After Heating');

% plot the temperature after cooling
subplot(2, 3, 4);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, T2);
h = colorbar;
xlabel(h, '[degC]');
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Temperature After Cooling');

% plot thermal dose
subplot(2, 3, 5);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, kdiff.cem43, [0, 1000]);
h = colorbar;
xlabel(h, '[CEM43]');
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Thermal Dose');

% plot lesion map
subplot(2, 3, 6);
imagesc(kgrid.y_vec * 1e3, kgrid.x_vec * 1e3, kdiff.lesion_map, [0, 1]);
colorbar;
ylabel('x-position [mm]');
xlabel('y-position [mm]');
axis image;
title('Ablated Tissue');

% set colormap and enlarge figure window
colormap(jet);