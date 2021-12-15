clear

%% ======================================================================== 
%%                        PARAMETERS 
%% ========================================================================

fc = 1.0e6;	           % fundamental frequency	[Hz]
p0 = 1.0e6;	           % initial pressure amplitude [Pa]
medium.c0 = 1500;      % speed of Sound [m/s]
medium.rho0 = 1000;    % density of medium [kg/m^3]

omegac = fc*2*pi;      % fundamental angular frequency
lambda = medium.c0/fc; % wavelength [m]

dx = lambda/16;        % step size in the x direction [m]
dt = 1/fc/8;           % timestep [s]

xlength = lambda*14;   % computational domain size in the x direction [m]
tlength = 60/fc;       % temporal domain size [s]
 
%% ======================================================================== 
%%                        COMPUTATIONAL DOMAIN
%% ========================================================================

mgrid = set_grid(dt, tlength, dx, xlength); 

%% ======================================================================== 
%%                        EXCITATION
%% ========================================================================
ts = [-4/fc:dt:4/fc].';
excit_ps = p0*sin(2*pi*fc*(ts)).*exp(-(ts).^2*fc^2/2);

%% ======================================================================== 
%%                        MEDIUM
%% ========================================================================
% 
speed_ratio = 1.0;
rho_ratio = 1.20;

medium.beta0 = 0;
medium.beta1 = 0;

medium.ca0 = 0.0;
medium.ca1 = 0.0;

z_1D = ([1:mgrid.num_x+1]-mgrid.num_x/2 - 1/2)*mgrid.dx;
medium.c    = medium.c0*(z_1D<=0) + medium.c0*speed_ratio*(z_1D>0);
medium.rho  = medium.rho0*(z_1D<=0) + medium.rho0*rho_ratio*(z_1D>0);
medium.ca   = medium.ca0*(z_1D<=0) + medium.ca1*(z_1D>0);
medium.beta = medium.beta0*(z_1D<=0) + medium.beta1*(z_1D>0);
% constant exponent b for the power law relation
medium.cb  = 2.0; 


%% ======================================================================== 
%%                        SIMULATION
%% ========================================================================
sensor_mask = ones(mgrid.num_x+1, 1);
reflection_order = 4;
tic
p_total = Forward1D(mgrid, medium, excit_ps, sensor_mask, reflection_order, 'correction');
toc
%% ======================================================================== 
%%                        RESULTS
%% ========================================================================

P_time(:) = p_total(:, end);

figure
plot(mgrid.t*1e6, P_time/1e6, 'b', 'LineWidth',2)

P_freq = abs(fftshift(fft(P_time)));
omega_ratio = mgrid.w/omegac;
P_dB = 20*log10(P_freq);
figure
plot(omega_ratio, P_dB, 'k', 'LineWidth',2)

p_end = P_time;
% save ('forward1D.mat', 'p_end', 'mgrid', 'medium')
