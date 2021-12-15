
clear
%% ======================================================================== 
%%                        PARAMETERS 
%% ========================================================================
fc = 1.0e6;	      % Fundamental frequency	[Hz]
p0 = 1;	      % Initial pressure amplitude  [Pa]
omegac = 2*fc*pi; % Fundamental angular frequency

medium.c0 = 1500;      % Speed of Sound [m/s]
medium.rho0 = 1000;    % Density of medium [kg/m^3]
medium.beta0 = 3.6;
medium.ca0 = 0.005;

lambda = medium.c0/fc; % Wave length [m]

dx = lambda/8;  %  increment in x axis   [m]
dy = lambda/8;  %  increment in z axis   [m]
dz = 14*lambda; %  increment in z axis   [m]

x_length = 60*lambda; % 30lambda
y_length = 60*lambda; % 30lambda
z_length = 14*lambda;

%% ======================================================================== 
%%                        COMPUTATIONAL DOMAIN
%% ========================================================================
mgrid = set_grid(0, 0, dx, x_length, dy, y_length, dz, z_length); 

%% ======================================================================== 
%%                        TRANSDUCER
%% ========================================================================

% TR_focus  = 8*lambda;    % transducer focus length [m]
% TR_radius = 5*lambda;    % Transducer diameter   [m]
% focus  = round(TR_focus/mgrid.dy);     % [grid points]
% num_el = round(TR_radius*2/mgrid.dx);    % [grid points]

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
imagesc(mgrid.y(round(mgrid.num_y/4):round(mgrid.num_y/4*3))*1e3,...
        mgrid.x(round(mgrid.num_x/4):round(mgrid.num_x/4*3))*1e3,...
        abs(source_p(round(mgrid.num_x/4):round(mgrid.num_x/4*3),...
            (round(mgrid.num_y/4):round(mgrid.num_y/4*3)))));
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
medium.NRL_gamma = 0.5;
medium.NRL_alpha = 0.8;
%% ======================================================================== 
%%                        SIMULATION
%% ========================================================================
% forward propagation of the fundamental pressure 
P_fundamental = Forward3D_fund(mgrid, medium, source_p, omegac, 0);

P_end = P_fundamental(:,:,end);
P_backward = Backward3D_fund(mgrid, medium, P_end, omegac, 0);
%% ======================================================================== 
%%                        RESULTS
%% ========================================================================

P_fundxz(:,:) = P_fundamental(:,round(mgrid.num_y/2), :);
P_fundyz(:,:) = P_fundamental(round(mgrid.num_x/2),:, :);
P_fundxy(:,:) = P_fundamental(:,:, end);

subplot(1,3,2)
imagesc(mgrid.y(round(mgrid.num_y/4):round(mgrid.num_y/4*3))*1e3,...
        mgrid.x(round(mgrid.num_x/4):round(mgrid.num_x/4*3))*1e3,...
        abs(P_fundxy((round(mgrid.num_x/4):round(mgrid.num_x/4*3)),...
                     (round(mgrid.num_y/4):round(mgrid.num_y/4*3))))...
        /max(max(abs(P_fundxy((round(mgrid.num_x/4):round(mgrid.num_x/4*3)),...
                             (round(mgrid.num_y/4):round(mgrid.num_y/4*3)))))));
xlabel ('y (mm)')
ylabel ('x (mm)')
colormap (jet)
axis image
colorbar


P_Bfundxz(:,:) = P_backward(:,round(mgrid.num_y/2), :);
P_Bfundyz(:,:) = P_backward(round(mgrid.num_x/2),:, :);
P_Bfundxy(:,:) = P_backward(:,:, 1);

subplot(1,3,3)
imagesc(mgrid.y(round(mgrid.num_y/4):round(mgrid.num_y/4*3))*1e3,...
        mgrid.x(round(mgrid.num_y/4):round(mgrid.num_y/4*3))*1e3,...
        abs(P_Bfundxy(round(mgrid.num_x/4):round(mgrid.num_x/4*3),...
                     (round(mgrid.num_y/4):round(mgrid.num_y/4*3))))...
        /max(max(abs(P_Bfundxy((round(mgrid.num_x/4):round(mgrid.num_x/4*3)),...
                               (round(mgrid.num_y/4):round(mgrid.num_y/4*3)))))));
xlabel ('y (mm)')
ylabel ('x (mm)')
colormap (jet)
axis image
colorbar