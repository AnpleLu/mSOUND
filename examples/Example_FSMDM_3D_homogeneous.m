
clear
% Simulation of a 3D homogeneous medium 
% using the frequency-specific mixed domain method

%% ======================================================================== 
%%                        PARAMETERS 
%% ========================================================================
fc = 1.0e6;	      % fundamental frequency	[Hz]
p0 = 1.0e6;	      % initial pressure amplitude  [Pa]
omegac = 2*fc*pi; % fundamental angular frequency

medium.c0 = 1500;      % speed of Sound [m/s]
medium.rho0 = 1000;    % density of medium [kg/m^3]
medium.beta0 = 3.6;    % nonlinearity coefficient 
medium.ca0 = 0.005;    %attenuation coefficient [dB/(MHz^y cm)] 

lambda = medium.c0/fc; % wavelength [m]

dx = lambda/8; % step size in the x direction [m]
dy = lambda/8; % step size in the y direction [m]
dz = lambda/8; % step size in the z direction [m]

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
TR_radius = 5*lambda;    % Transducer radius [m]
focus  = round(TR_focus/mgrid.dz);     % [grid points]
num_el = round(TR_radius*2/mgrid.dx);  % [grid points]

%% ======================================================================== 
%%                        EXCITATION
%% ========================================================================
X = ones(mgrid.num_y,1)*mgrid.x;
Y = (ones(mgrid.num_x,1)*mgrid.y).';

% calculate the delay for each element [m]   
delay = sqrt(X.^2 + Y.^2 + TR_focus^2)/medium.c0;   
delay = delay - min(min(delay));

source_p = p0*exp(1i*omegac*delay);

RHO = sqrt(X.^2+Y.^2);
source_p(RHO>TR_radius) = 0;
clear X Y RHO

imagesc(abs(source_p))
%% ======================================================================== 
%%                        MEDIUM
%% ========================================================================    
medium.c    = medium.c0;
medium.rho  = medium.rho0;
medium.beta = medium.beta0;
medium.ca   = medium.ca0;
medium.cb   = 2.0;  %power law exponent

%% ======================================================================== 
%%                        SIMULATION
%% ========================================================================
% forward propagation of the fundamental pressure 
P_fundamental = Forward3D_fund(mgrid, medium, source_p, omegac, 0);

% forward propagation of the second-harmonic pressure   
P_second      = Forward3D_sec(mgrid, medium, P_fundamental, omegac);

%% ======================================================================== 
%%                        RESULTS
%% ========================================================================

P_fundxz(:,:) = P_fundamental(40:120,round(mgrid.num_y/2), :);
P_fundyz(:,:) = P_fundamental(round(mgrid.num_x/2), 40:120, :);
P_fundxy(:,:) = P_fundamental(40:120,40:120, focus+1);
P_secxz(:,:)  = P_second(40:120,round(mgrid.num_y/2), :);
P_secyz(:,:)  = P_second(round(mgrid.num_x/2), 40:120, :);
P_secxy(:,:)  = P_second(40:120,40:120, focus+1);
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
ylabel ('x(mm)')
colormap (jet)
axis image
colorbar
subplot(1,3,3)
imagesc(mgrid.y(40:120)*1e3, mgrid.x(40:120)*1e3,abs(P_fundxy)/max(max(abs(P_fundxy))));
xlabel ('z(mm)')
ylabel ('x(mm)')
colormap (jet)
axis image
colorbar

figure
subplot(1,3,1)
imagesc(mgrid.z*1e3, mgrid.x(40:120)*1e3,abs(P_secxz)/max(max(abs(P_secxz))));
xlabel ('z(mm)')
ylabel ('x(mm)')
colormap (jet) 
axis image
colorbar
subplot(1,3,2)
imagesc(mgrid.z*1e3, mgrid.x(40:120)*1e3,abs(P_secyz)/max(max(abs(P_secyz))));
xlabel ('z(mm)')
ylabel ('x(mm)')
colormap (jet)
axis image
colorbar
subplot(1,3,3)
imagesc(mgrid.y(40:120)*1e3, mgrid.x(40:120)*1e3,abs(P_secxy)/max(max(abs(P_secxy))));
xlabel ('z(mm)')
ylabel ('x(mm)')
colormap (jet)
axis image
colorbar





