clear
% Reducing the spatial aliasing error using the non-reflecting layer
%% ======================================================================== 
%%                        PARAMETERS 
%% ========================================================================
fc = 1.0e6;	           % fundamental frequency	[Hz]
p0 = 0.6e6;	           % initial pressure amplitude  [Pa]
medium.c0 = 1500;      % speed of Sound [m/s]
medium.rho0 = 1000;    % density of medium [kg/m^3]

omegac = fc*2*pi;      % fundamental angular frequency
lambda = medium.c0/fc; % wavelength [m]

dx = lambda/6;         % step size in the x direction [m]
dy = lambda/6;         % step size in the y direction [m]

x_length = lambda*25+dx; % computational domain size in the x direction [m]
y_length = lambda*20;    % computational domain size in the y direction [m]

%% ======================================================================== 
%%                        COMPUTATIONAL DOMAIN
%% ========================================================================

mgrid = set_grid(0,0, dx, x_length, dy, y_length);  

%% ======================================================================== 
%%                        TRANSDUCER
%% ========================================================================

TR_focus  = 12*lambda;  % transducer focal length [m]
TR_radius = 6*lambda;   % Transducer radius [m]

focus  = round(TR_focus/mgrid.dy);   % [grid points]

%% ======================================================================== 
%%                        EXCITATION
%% ========================================================================

% calculate the phase delay
delay = sqrt((mgrid.x).^2 + (TR_focus)^2)/medium.c0;  % [s]
delay = delay - min(delay);

% define the excitation pulse 
excit_p = p0*exp(1i*omegac*delay);  
excit_p(abs(mgrid.x)>TR_radius) = 0;
%% ======================================================================== 
%%                        MEDIUM
%% ========================================================================
medium.c    = 1500;   % speed of sound [m/s]      
medium.rho  = 1000;   % density [kg/m^3]             
medium.beta = 3.6;    % nonlinear coefficient     
medium.ca   = 0;      % attenuation coefficient [dB/(MHz^y cm)]     
medium.cb   = 2.0;    % power law exponent    

medium.NRL_gamma = 0.8;  % 0.1; 0.5; 0.8
medium.NRL_alpha = 0.05; % 0.5; 0.2; 0.05

%% ======================================================================== 
%%                        SIMULATION
%% ========================================================================

% forward propagation of the fundamental pressure 
P_fundamental = Forward2D_fund(mgrid, medium, excit_p, omegac, 0, 'NRL');

% forward propagation of the second-harmonic pressure     
P_second = Forward2D_sec(mgrid, medium, P_fundamental, omegac);
%% ======================================================================== 
%%                        RESULTS
%% ========================================================================

P_fund = P_fundamental(:,2:end);
P_secd = P_second(:,2:end);

figure
py = P_fund(round(mgrid.num_x/2), :);
plot(mgrid.y(1:end)*1e3 - min(mgrid.y(1:end)*1e3), abs(py)/1e6, 'k','LineWidth',2)
