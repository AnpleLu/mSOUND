
clear
% One-step reconstruction of the source pressure 
% distribution with FSMDM in a 3D homogeneous medium
%% ======================================================================== 
%%                        PARAMETERS 
%% ========================================================================
fc = 1.0e6;	      % fundamental frequency	[Hz]
p0 = 1;	          % initial pressure amplitude  [Pa]
omegac = 2*fc*pi; % fundamental angular frequency

medium.c0 = 1500;      % speed of Sound [m/s]
medium.rho0 = 1000;    % density of medium [kg/m^3]
medium.beta0 = 3.6;
medium.ca0 = 0.005;

lambda = medium.c0/fc; % wavelength [m]

dx = lambda/8;  % step size in the x direction [m]
dy = lambda/8;  % step size in the y direction [m]
dz = 0.021;  % step size in the z direction [m]

x_length = 30*lambda;  % computational domain size in the x direction [m]
y_length = 30*lambda;  % computational domain size in the y direction [m]
z_length = 0.021;  % computational domain size in the z direction [m]

%% ======================================================================== 
%%                        COMPUTATIONAL DOMAIN
%% ========================================================================

mgrid = set_grid(0, 0, dx, x_length, dy, y_length, dz, z_length); 

%% ======================================================================== 
%%                        TRANSDUCER
%% ========================================================================

TR_focus  = 8*lambda;    % transducer focal length [m]
TR_radius = 5*lambda;    % transducer radius [m]
focus  = round(TR_focus/mgrid.dy);     % [grid points]
num_el = round(TR_radius*2/mgrid.dx);  % [grid points]

%% ======================================================================== 
%%                        EXCITATION
%% ========================================================================
X = ones(mgrid.num_y,1)*mgrid.x;
Y = (ones(mgrid.num_x,1)*mgrid.y).';
RHO = sqrt(X.^2+Y.^2);
RHO(RHO>2.5*lambda) = 0;

source_map1 = sqrt((X-lambda*3).^2 + (Y+lambda*3).^2);   
source_map2 = sqrt((X+lambda*3).^2 + (Y+lambda*3).^2); 
source_map3 = sqrt((X-lambda*3).^2 + (Y-lambda*3).^2); 
source_map4 = sqrt((X+lambda*3).^2 + (Y-lambda*3).^2); 
TR_radius = 2.5*lambda; 
source_map1(source_map1<=TR_radius) = 1;
source_map2(source_map2<=TR_radius) = 1;
source_map3(source_map3<=TR_radius) = 1;
source_map4(source_map4<=TR_radius) = 1;
source_map = source_map1  + source_map2 + source_map3 + source_map4;

source_map(source_map<1)=0;
source_p = source_map;

clear excit_ptemp X Y

figure
subplot(1,3,1)
imagesc(mgrid.y*1e3, mgrid.x*1e3,abs(source_p));
xlabel ('y (mm)')
ylabel ('x (mm)')
colormap (jet)
axis image
colorbar
%% ======================================================================== 
%%                        MEDIUM
%% ========================================================================    
medium.c    = medium.c0;
medium.rho  = medium.rho0;
medium.beta = medium.beta0;
medium.ca   = medium.ca0;
medium.cb   = 2.0;
medium.NRL_gamma = 0.1;
medium.NRL_alpha = 0.1;
%% ======================================================================== 
%%                        SIMULATION
%% ========================================================================
% forward propagation of the fundamental pressure 
P_fundamental = Forward3D_fund(mgrid, medium, source_p, omegac, 0, 'NRL');

P_end = P_fundamental(:,:,end);
P_backward = Backward3D_fund(mgrid, medium, P_end, omegac, 0, 'NRL');
%% ======================================================================== 
%%                        RESULTS
%% ========================================================================

P_fundxz(:,:) = P_fundamental(:,round(mgrid.num_y/2), :);
P_fundyz(:,:) = P_fundamental(round(mgrid.num_x/2),:, :);
P_fundxy(:,:) = P_fundamental(:,:, end);

subplot(1,3,2)
imagesc(mgrid.y*1e3, mgrid.x*1e3,abs(P_fundxy)/max(max(abs(P_fundxy))));
xlabel ('y (mm)')
ylabel ('x (mm)')
colormap (jet)
axis image
colorbar

P_Bfundxz(:,:) = P_backward(:,round(mgrid.num_y/2), :);
P_Bfundyz(:,:) = P_backward(round(mgrid.num_x/2),:, :);
P_Bfundxy(:,:) = P_backward(:,:, 1);

subplot(1,3,3)
imagesc(mgrid.y*1e3, mgrid.x*1e3,abs(P_Bfundxy)/max(max(abs(P_Bfundxy))));
xlabel ('y (mm)')
ylabel ('x (mm)')
colormap (jet)
axis image
colorbar





