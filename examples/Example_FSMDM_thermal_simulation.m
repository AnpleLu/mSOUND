clear
% Integrating mSOUND with k-Wave for thermal simulations

%% ======================================================================== 
%%                        PARAMETERS 
%% ========================================================================
fc = 1.0e6;	           % Fundamental frequency	[Hz]
p0 = 0.8e6;	           % Initial pressure amplitude  [Pa]
medium.c0 = 1500;      % Speed of Sound [m/s]

omegac = fc*2*pi;      % Fundamental angular frequency
lambda = medium.c0/fc; % Wave length [m]

dx = lambda/10;
dy = lambda/10;

x_length = lambda*30+dx;
y_length = lambda*30;

%% ======================================================================== 
%%                        COMPUTATIONAL DOMAIN
%% ========================================================================

mgrid = set_grid(0, 0, dx, x_length, dy, y_length);  %% wavevector and x, y

%% ======================================================================== 
%%                        TRANSDUCER
%% ========================================================================

TR_focus  = 24*lambda;  % transducer focus length [m]
TR_radius = 6*lambda;   % Transducer diameter   [m]

%% ======================================================================== 
%%                        EXCITATION
%% ========================================================================

% calculate the phase delay
delay = sqrt((mgrid.x).^2 + (TR_focus)^2)/medium.c0;  % [s]
delay = delay - min(delay);

source_p = p0*exp(1i*omegac*delay);  % define the pulse 
source_p(abs(mgrid.x)>TR_radius) = 0;
%% ======================================================================== 
%%                        MEDIUM
%% ========================================================================
medium.c    = 1500;   % speed of sound [m/s]      
medium.rho  = 1000;   % density [kg/m^3]             
medium.beta = 0;      % nonlinear coefficient     
medium.ca   = 0.75;   % attenuation coefficient [dB/(MHz^y cm)]     
medium.cb   = 1.5;    % power law exponent    

medium.NRL_gamma = 0.5;
medium.NRL_alpha = 0.05;
%% ======================================================================== 
%%                        SIMULATION
%% ========================================================================
% forward propagation of the fundamental pressure 
P_fundamental = Forward2D_fund(mgrid, medium, source_p, omegac, 0, 'NRL');
p = abs(P_fundamental(:,2:end));

% figure
% plot(mgrid.y*1e3, p(round(mgrid.num_x/2), :), 'b', 'LineWidth', 2);
%% ======================================================================== 
%%                        thermal simulation
%% ========================================================================

% convert the absorption coefficient to nepers/m
alpha_np = db2neper(medium.ca, medium.cb) * (omegac).^medium.cb;

Q = alpha_np .* p.^2 ./ (medium.rho .* medium.c);
% set the background temperature and heating term
source.Q = Q;
source.T0 = 37;

clear medium
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
% scaleFig(1.5, 1);