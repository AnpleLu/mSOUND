clear
%% ======================================================================== 
%%                        PARAMETERS 
%% ========================================================================
fc = 1.0e6;	           % fundamental frequency	[Hz]
p0 = 2.0e6;	           % initial pressure amplitude  [Pa]
medium.c0 = 1479;      % speed of Sound [m/s]
medium.rho0 = 1000;    % density of medium [kg/m^3]

omegac = fc*2*pi;      % fundamental angular frequency
lambda = medium.c0/fc; % wavelength [m]

dx = lambda/8;         % step size in the x direction [m]
dy = lambda/8;         % step size in the y direction [m]

x_length = lambda*50;  % computational domain size in the x direction [m]
y_length = lambda*22;  % computational domain size in the y direction [m]

%% ======================================================================== 
%%                        COMPUTATIONAL DOMAIN
%% ========================================================================

mgrid = set_grid(0,0, dx, x_length, dy, y_length); 

%% ======================================================================== 
%%                        TRANSDUCER
%% ========================================================================

TR_focus  = 14*lambda;     % transducer focal length [m]
TR_radius = 7.5*lambda;    % transducer radius [m]

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
load tissue_model.mat
medium.c    = c;
medium.rho  = rho;
medium.beta = beta;
medium.ca   = ca;      % power law relation constant     
medium.cb   = cb; 

medium.NRL_gamma = 0.5;
medium.NRL_alpha = 0.1;
%% ======================================================================== 
%%                        SIMULATION
%% ========================================================================

% the maximum order of reflection included in the simulation  
reflection_order = 2;

% forward propagation of the fundamental pressure 
P_fundamental = Forward2D_fund(mgrid, medium, excit_p, omegac, reflection_order, 'NRL');

% forward propagation of the second-harmonic pressure     
P_second = Forward2D_sec(mgrid, medium, P_fundamental, omegac, 'NRL');

% forward propagation of the second-harmonic pressure     
P_second2 = Forward2D_sec(mgrid, medium, P_fundamental, omegac, 'NRL', 'correction');
%% ======================================================================== 
%%                        RESULTS
%% ========================================================================

P_fund = P_fundamental(:,2:end);
P_secd = P_second(:,2:end);
P_secd2 = P_second2(:,2:end);
figure
subplot(1,2,1)
imagesc(mgrid.y*1000,-mgrid.x(100:300)*1000,abs(P_fund(100:300, :)));
xlabel ('y(mm)')
ylabel ('x(mm)')
colormap (jet) 
axis image
% colorbar
subplot(1,2,2)
imagesc(mgrid.y*1000,-mgrid.x(100:300)*1000,abs(P_secd(100:300,:)));
xlabel ('y(mm)')
ylabel ('x(mm)')
colormap (jet)
axis image
% colorbar

figure
imagesc(mgrid.y*1000,-mgrid.x(100:300)*1000,abs(P_secd2(100:300,:)));
xlabel ('y(mm)')
ylabel ('x(mm)')
colormap (jet)
axis image
% colorbar


