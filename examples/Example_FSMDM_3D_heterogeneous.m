
clear
% Simulation of a 3D heterogeneous medium 
% using the frequency-specific mixed domain method

%% ======================================================================== 
%%                        PARAMETERS 
%% ========================================================================
fc = 1.0e6;	      % fundamental frequency	[Hz]
p0 = 1.0e6;	      % initial pressure amplitude  [Pa]
omegac = 2*fc*pi; % fundamental angular frequency

medium.c0 = 1500;      % speed of Sound [m/s]
medium.rho0 = 1000;    % density of medium [kg/m^3]
medium.beta0 = 3.6;
medium.ca0 = 0.005;

lambda = medium.c0/fc; % wavelength [m]

dx = lambda/8;  % step size in the x direction [m]
dy = lambda/8;  % step size in the y direction [m]
dz = lambda/8;  % step size in the z direction [m]

x_length = 20*lambda;  % computational domain size in the x direction [m]
y_length = 20*lambda;  % computational domain size in the y direction [m]
z_length = 14*lambda;  % computational domain size in the z direction [m]
%% ======================================================================== 
%%                        COMPUTATIONAL DOMAIN
%% ========================================================================

mgrid = set_grid(0, 0, dx, x_length, dy, y_length, dz, z_length); 

%% ======================================================================== 
%%                        TRANSDUCER
%% ========================================================================

TR_focus  = 8*lambda;    % transducer focal length [m]
TR_radius = 5*lambda;    % transducer radius  [m]
focus  = round(TR_focus/mgrid.dy);     % [grid points]
num_el = round(TR_radius*2/mgrid.dx);  % [grid points]

%% ======================================================================== 
%%                        EXCITATION
%% ========================================================================
X = ones(mgrid.num_y,1)*mgrid.x;
Y = (ones(mgrid.num_x,1)*mgrid.y).';
delay = sqrt(X.^2 + Y.^2 + TR_focus^2)/medium.c0;   % calculate the delay for each element [m]   
delay = delay - min(min(delay));

source_p = p0*exp(1i*omegac*delay);

RHO = sqrt(X.^2+Y.^2);
source_p(RHO>TR_radius) = 0;
clear excit_ptemp X Y

%% ======================================================================== 
%%                        MEDIUM
%% ========================================================================    
speed_water = 1500*1.00;
speed_fat   = 1500*1.50;

rho_water = 1000*1.00;
rho_fat   = 1000*1.50;
 
cb       = 2.0;      % constant exponent b for the power law relation
ca_water = 0.005;
ca_fat   = 0.005;

beta_water = 3.6;
beta_fat   = 3.6;
 
[c1, rho1, beta1, ca1] = skull3D(mgrid.dx, mgrid.dy, mgrid.dz, ...
                         mgrid.num_x, mgrid.num_y,mgrid.num_z, ...
                         fc, medium.c0,...
                         speed_water, speed_fat, ...
                         rho_water, rho_fat,...
                         beta_water, beta_fat, ...
                         ca_water, ca_fat, 0);
                     
medium.c    = medium.c0*ones(mgrid.num_x, mgrid.num_y, mgrid.num_z+1);
medium.rho  = medium.rho0*ones(mgrid.num_x, mgrid.num_y, mgrid.num_z+1);
medium.beta = medium.beta0*ones(mgrid.num_x, mgrid.num_y, mgrid.num_z+1);
medium.ca   = medium.ca0*ones(mgrid.num_x, mgrid.num_y, mgrid.num_z+1);

medium.c(:,:,2:end)    = c1;
medium.rho(:,:,2:end)  = rho1;
medium.beta(:,:,2:end) = beta1;
medium.ca(:,:,2:end)   = ca1;       % power law relation constant     
medium.cb = 2.0;

c = medium.c;
rho = medium.rho;
beta = medium.beta;
ca = medium.ca;
save ('3D_strongly_heterogeneous_model.mat', 'c', 'rho','beta', 'ca')

% load 3D_strong_heterogeneous_model.mat
% medium.c = c;
% medium.rho = rho;
% medium.beta = beta;
% medium.ca = ca;
% medium.cb = 2.0;

medium.NRL_gamma = 0.5;
medium.NRL_alpha = 0.1;
figure
imagesc(mgrid.z*1e3, mgrid.y*1e3, squeeze(medium.c(round(mgrid.num_x/2),:,:)))
axis image
colormap (jet)
%% ======================================================================== 
%%                        SIMULATION
%% ========================================================================
reflection_order = 4;
P_fundamental = Forward3D_fund(mgrid, medium, source_p, omegac,...
                reflection_order, 'NRL', 'correction');
% P_second      = Forward3D_sec(mgrid, medium, P_fundamental, omegac,...
%                 'NRL', 'correction');

%% ======================================================================== 
%%                        RESULTS
%% ========================================================================

P_fundxz(:,:) = P_fundamental(40:120,round(mgrid.num_y/2), :);
P_fundyz(:,:) = P_fundamental(round(mgrid.num_x/2), 40:120, :);
P_fundxy(:,:) = P_fundamental(40:120,40:120, focus+1);
% P_secxz(:,:)  = P_second(40:120,round(mgrid.num_y/2), :);
% P_secyz(:,:)  = P_second(round(mgrid.num_x/2), 40:120, :);
% P_secxy(:,:)  = P_second(40:120,40:120, focus+1);
figure
subplot(1,3,1)
imagesc(mgrid.z*1e3, mgrid.x(40:120)*1e3, abs(P_fundxz)/max(max(abs(P_fundxz))));
xlabel ('z(mm)')
ylabel ('x(mm)')
colormap (jet) 
axis image
colorbar
subplot(1,3,2)
imagesc(mgrid.z*1e3, mgrid.y(40:120)*1e3,abs(P_fundyz)/max(max(abs(P_fundyz))));
xlabel ('z(mm)')
ylabel ('y(mm)')
colormap (jet)
axis image
colorbar
subplot(1,3,3)
imagesc(mgrid.y(40:120)*1e3, mgrid.x(40:120)*1e3,abs(P_fundxy)/max(max(abs(P_fundxy))));
xlabel ('y(mm)')
ylabel ('x(mm)')
colormap (jet)
axis image
colorbar

% subplot(2,3,4)
% imagesc(mgrid.z*1e3, mgrid.x(40:120)*1e3,abs(P_secxz)/max(max(abs(P_secxz))));
% xlabel ('z(mm)')
% ylabel ('x(mm)')
% colormap (jet) 
% axis image
% colorbar
% subplot(2,3,5)
% imagesc(mgrid.z*1e3, mgrid.x(40:120)*1e3,abs(P_secyz)/max(max(abs(P_secyz))));
% xlabel ('z(mm)')
% ylabel ('x(mm)')
% colormap (jet)
% axis image
% colorbar
% subplot(2,3,6)
% imagesc(mgrid.y(40:120)*1e3, mgrid.x(40:120)*1e3,abs(P_secxy)/max(max(abs(P_secxy))));
% xlabel ('z(mm)')
% ylabel ('x(mm)')
% colormap (jet)
% axis image
% colorbar





