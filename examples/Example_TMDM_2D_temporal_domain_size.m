clear
% Selecting the proper temporal domain size for the TMDM

%% ======================================================================== 
%%                        PARAMETERS 
%% ========================================================================
fc = 1.0e6;	           % fundamental frequency	[Hz]
p0 = 1.0e6;	           % initial pressure amplitude  [Pa]
medium.c0 = 1500;      % speed of Sound [m/s]
medium.rho0 = 1000;    % density of medium [kg/m^3]

omegac = fc*2*pi;      % fundamental angular frequency
lambda = medium.c0/fc; % wavelength [m]

dx = lambda/8;         % step size in the x direction [m]
dy = lambda/8;         % step size in the y direction [m]
dt = 1/fc/16;          % timestep  [s]

x_length = lambda*50;  % computational domain size in the x direction [m]
y_length = lambda*22;  % computational domain size in the y direction [m]
% t_length = 4/fc*3.5; % 5.4   4.5   3.5
t_length = 3.2e-5; % 8e-6
%% ======================================================================== 
%%                        COMPUTATIONAL DOMAIN
%% ========================================================================

mgrid = set_grid(dt, t_length, dx, x_length, dy, y_length);  

%% ======================================================================== 
%%                        TRANSDUCER
%% ========================================================================

TR_focus  = 14*lambda;     % transducer focal length [m]
TR_radius = 7.5*lambda;    % Transducer radius [m]
focus  = round(TR_focus/mgrid.dy);  % [grid points]

% calculate the phase delay
delay = sqrt((mgrid.x).^2 + (TR_focus)^2)/medium.c0;  % [s]
delay = delay - min(delay);

%% ======================================================================== 
%%                        EXCITATION
%% ========================================================================

% define the pulse length
ts = [-4/fc:dt:4/fc].';
delay = repmat(delay, length(ts),1);
ts = repmat(ts, 1,mgrid.num_x);

% define the excitation pulse
source_p = p0*sin(2*pi*fc*(ts+delay)).*exp(-(ts+delay).^2*fc^2/2);
source_p(:, abs(mgrid.x)>TR_radius) = 0;
%% ======================================================================== 
%%                        MEDIUM
%% ========================================================================
load tissue_model.mat
medium.c    = c;
medium.rho  = rho;
medium.beta = beta;
medium.ca   = ca;       % power law relation constant     
medium.cb   = cb;

% visualize the 2D model
figure
imagesc(mgrid.y*1e3, mgrid.x*1e3,beta)
colormap(jet)
axis image
hold on
plot((TR_focus+min(mgrid.y))*1000,0,'r.','MarkerSize',20)
graph1 = line([min(mgrid.y)*1e3,min(mgrid.y)*1e3],[-TR_radius*1000,TR_radius*1000], 'Color','red');
set(graph1,'LineWidth',4);
axis equal
%% ======================================================================== 
%%                        SIMULATION
%% ========================================================================
% define where the pressure will be recorded
sensor_mask = zeros(mgrid.num_x, mgrid.num_y+1);
sensor_mask(round(mgrid.num_x/2), 113) = 1;

% the maximum order of reflection included in the simulation
reflection_order = 2; 

% 2D forward simulation
tic
p_time = Forward2D(mgrid, medium, source_p, sensor_mask,...
         reflection_order, 'animation');
toc
% figure
% subplot (2,2,2)
% plot(mgrid.t*1e6, p_time/1e6, 'k', 'LineWidth', 2)
% xlabel ('Time (\mus)')
% ylabel ('Pressure (MPa)')
% 
% figure
% plot(mgrid.w/2/pi/1e6,20*log10(abs(fftshift(fft(p_time)))), 'k', 'LineWidth', 2)
% save ('animation2D_c.mat', 'mgrid', 'p_time')
% 
