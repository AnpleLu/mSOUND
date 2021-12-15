clear
% Shock wave simulations with TMDM
%% ======================================================================== 
%%                        PARAMETERS 
%% ========================================================================
fc = 1.0e6;	           % fundamental frequency	[Hz]
p0 = 6.5e6;	           % initial pressure amplitude  [Pa]
medium.c0 = 1500;      % speed of Sound [m/s]
medium.rho0 = 1000;    % density of medium [kg/m^3]

omegac = fc*2*pi;      % fundamental angular frequency
lambda = medium.c0/fc; % wavelength [m]

dx = lambda/8;     % step size in the x direction [m]
dy = lambda/200;   % step size in the y direction [m], 1000
dt = 1/fc/200;     % timestep [s] 1000

x_length = lambda*30+dx;
y_length = lambda*12;
t_length = 16/fc; 

%% ======================================================================== 
%%                        COMPUTATIONAL DOMAIN
%% ========================================================================

mgrid = set_grid(dt, t_length, dx, x_length, dy, y_length);  %% wavevector and x, y

%% ======================================================================== 
%%                        TRANSDUCER
%% ========================================================================

TR_focus  = 12*lambda;  % transducer focal length [m]
TR_radius = 6*lambda;   % transducer radius [m]
focus  = round(TR_focus/mgrid.dy);   % [grid points]

% calculate the phase delay
delay = sqrt((mgrid.x).^2 + (TR_focus)^2)/medium.c0;  % [s]
delay = delay - min(delay);

%% ======================================================================== 
%%                        EXCITATION
%% ========================================================================
% define the pulse length
ts = [0:dt:16/fc].';
delay = repmat(delay, length(ts),1);
ts = repmat(ts, 1,mgrid.num_x);

% define the excitation pulse
source_p = p0*sin(2*pi*fc*(ts+delay));
source_p(:, abs(mgrid.x)>TR_radius) = 0;

%% ======================================================================== 
%%                        MEDIUM
%% ========================================================================

medium.c    = 1500;
medium.rho  = 1000;
medium.beta = 3.6;
medium.ca   = 0.02; % 0.1  % power law relation constant     
medium.cb   = 2.0;

medium.NRL_gamma = 0.5;
medium.NRL_alpha = 0.05;
%% ======================================================================== 
%%                        SIMULATION
%% ========================================================================
% define where the pressure will be recorded
sensor_mask = zeros(mgrid.num_x, mgrid.num_y+1);
sensor_mask(round(mgrid.num_x/2), focus+1) = 1;

% 2D forward simulation
tic
p_time = Forward2D(mgrid, medium, source_p, sensor_mask, 0);
toc
%% ======================================================================== 
%%                        RESULTS
%% ========================================================================

p_focus = squeeze(p_time);

figure
plot(mgrid.t*1e6, p_focus/1e6, 'r', 'LineWidth',2)
xlabel ('Time (\mus)')
ylabel ('Pressure (MPa)')

figure
plot(mgrid.w/2/pi/1e6,20*log10(abs(fftshift(fft(p_focus)))), 'r', 'LineWidth',2)
xlabel ('Frequency (MHz)')
ylabel ('Pressure (dB)')

% save ('shock_smallatt.mat', 'p_focus', 'mgrid')
