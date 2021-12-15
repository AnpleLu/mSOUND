clear
% Simulation of a 2D homogeneous medium 
% using the frequency-specific mixed domain method 
%% ======================================================================== 
%%                        PARAMETERS 
%% ========================================================================
fc = 1.0e6;	           % fundamental frequency	[Hz]
p0 = 0.6e6;	           % initial pressure amplitude  [Pa]
medium.c0 = 1500;      % speed of Sound [m/s]
medium.rho0 = 1000;    % density of medium [kg/m^3]

omegac = fc*2*pi;      % fundamental angular frequency
lambda = medium.c0/fc; % wavelength [m]

dx = lambda/6;    % step size in the x direction [m]
dy = lambda/6;    % step size in the y direction [m]

x_length = lambda*25+dx; % computational domain size in the x direction [m]
y_length = lambda*20;    % computational domain size in the y direction [m]

%% ======================================================================== 
%%                        COMPUTATIONAL DOMAIN
%% ========================================================================

mgrid = set_grid(0, 0, dx, x_length, dy, y_length);

%% ======================================================================== 
%%                        TRANSDUCER
%% ========================================================================

TR_focus  = 12*lambda;  % transducer focal length [m]
TR_radius = 6*lambda;   % transducer radius [m]

focus  = round(TR_focus/mgrid.dy);   % [grid points]

%% ======================================================================== 
%%                        EXCITATION
%% ========================================================================

% calculate the phase delay
delay = sqrt((mgrid.x).^2 + (TR_focus)^2)/medium.c0;  % [s]
delay = delay - min(delay);

excit_p = p0*exp(1i*omegac*delay);  % define the pulse 
excit_p(abs(mgrid.x)>TR_radius) = 0;
%% ======================================================================== 
%%                        MEDIUM
%% ========================================================================
medium.c    = 1500;   % speed of sound [m/s]      
medium.rho  = 1000;   % density [kg/m^3]             
medium.beta = 3.6;    % nonlinear coefficient     
medium.ca   = 0;      % attenuation coefficient [dB/(MHz^y cm)]     
medium.cb   = 2.0;    % power law exponent    


medium.NRL_gamma = 0.5; % constant for non-reflecting layer 
medium.NRL_alpha = 0.1; % decay factor for non-relfecting layer 
%% ======================================================================== 
%%                        SIMULATION
%% ========================================================================
% forward propagation of the fundamental pressure 
P_fundamental = Forward2D_fund(mgrid, medium, excit_p, omegac, 0, 'NRL');

% forward propagation of the second-harmonic pressure     
P_second = Forward2D_sec(mgrid, medium, P_fundamental, omegac, 'NRL');
%% ======================================================================== 
%%                        RESULTS
%% ========================================================================

P_fund = P_fundamental(:,2:end);
P_secd = P_second(:,2:end);

figure
subplot(1,2,1)
imagesc(mgrid.y*1000,mgrid.x*1000,abs(P_fund)/1e6);
xlabel ('y(mm)')
ylabel ('x(mm)')
colormap (jet) 
axis image
h = colorbar;
xlabel(h, '[MPa]');
subplot(1,2,2)
imagesc(mgrid.y*1000,mgrid.x*1000,abs(P_secd)/1e6);
xlabel ('y(mm)')
ylabel ('x(mm)')
colormap (jet)
axis image
h = colorbar;
xlabel(h, '[MPa]');

