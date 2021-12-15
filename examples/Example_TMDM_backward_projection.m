
clear
% Photoacoustic image reconstruction using backward projection

%% ======================================================================== 
%%                        k-Wave SIMULATION
%% ========================================================================
Nx = 800;
Ny = 300;

fc = 0.7e6;      % fundamental frequency [Hz]
c0 = 1500;       % speed of Sound [m/s]
p0 = 1;          % pressure magnitude [Pa]
lambda = c0/fc;  % wavelength [m]

dx = 1.3393e-04;  % step size in the x direction [m]
dy = 1.3393e-04;  % step size in the y direction [m]

kgrid = makeGrid(Nx, dx, Ny, dy);
%% create the time array
dt = 1.1161e-08;
cfl = c0*dt/dx;
display(cfl);

kgrid.t_array = 0:dt:3.4286e-05; 

load mSOUND_logo.mat
source.p = p0;
%% assign the properties of the propagation medium
medium.sound_speed = 1500*ones(Nx, Ny);
medium.BonA        = 2.0*(0-1.0);
medium.density     = 1000*ones(Nx, Ny);
medium.alpha_coeff = 0.0;   % power law absorption coefficient  
medium.alpha_power = 2.0;   % power law exponent

% c_ill = medium.sound_speed;
% c_ill(s==1) = 1600;
% figure
% subplot(1,3,1)
% imagesc(-kgrid.y_vec*1e3 + max(kgrid.y_vec*1e3), kgrid.x_vec*1e3, fliplr(c_ill))
% axis image
% colormap(jet);
% hold on
% graph1 = line([-min(kgrid.y_vec)*1e3 + round(Ny_p/2)*dy*1e3,...
%                -min(kgrid.y_vec)*1e3 + round(Ny_p/2)*dy*1e3],...
%               [min(kgrid.x_vec)*1000,max(kgrid.x_vec)*1000], 'Color','red');
% set(graph1,'LineWidth',2);
% axis equal
% hold off

%% define source and sensor
source.p_mask = s;
sensor.mask   = zeros(Nx, Ny);
offcenter     = floor(200/2); % propagation direction

sensor.mask(:, 250) = 1; 

%% excitation signal
sample_freq = 1/dt; % [Hz]
signal_freq = 1e6;  % [Hz]
source_mag = p0;    % [Pa] 

%% set the input options
input_args    = {'DisplayMask', sensor.mask, 'DataCast', 'single',...
                 'PlotScale',[-source_mag/2, source_mag/2]};
sensor.record = {'p'};

% run the simulation
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

timelength = length(kgrid.t_array);
sensor_data.p = reshape(sensor_data.p, Nx, timelength);
kwave_P_time = sensor_data.p;

% reset the initial pressure
source.p0 = 0;
source.p = 0;

% assign the time reversal data
sensor.time_reversal_boundary_data = sensor_data.p;
clear sensor_data 
%% ======================================================================== 
%%                        TIME REVERSAL
%% ========================================================================
% run the time reversal reconstruction
input_args    = {'DisplayMask', sensor.mask,'DataCast', 'single',...
                 'PlotScale',[-source_mag/2, source_mag/2]};
p0_recon = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

%% frequency domain
figure

subplot(1,3,2)
imagesc(kgrid.y_vec(kgrid.Ny/2-offcenter+1:kgrid.Ny/2 + offcenter)*1e3 -...
        min(kgrid.y_vec(kgrid.Ny/2-offcenter+1:kgrid.Ny/2 + offcenter)*1e3 ), kgrid.x_vec*1e3,...
       (p0_recon(:,kgrid.Ny/2-offcenter+1:kgrid.Ny/2 + offcenter)));
colormap(jet);
axis image;
colorbar;
caxis ([-0.05, 0.13])

%% ======================================================================== 
%%                        BACKWARD PROJECTION
%% ========================================================================
medium.c0 = 1500;   % speed of Sound [m/s]
dx = lambda/16;     % step size in the x direction [m]
dy = lambda/16;     % step size in the x direction [m]
num_x = 800;
num_y = 200;

omegac = fc*2*pi; % fundamental angular frequency

dt = 1/fc/16;     % timestep [s]

x_length = num_x*dx;
y_length = num_y*dy;
t_length = 24/fc; 

%% ============================================================
%%                   COMPUTATIONAL DOMAIN
%% ============================================================

mgrid = set_grid(dt, t_length, dx, x_length, dy, y_length);  %% wavevector and x, y

%% ==========================================================
%%                   EXCITATION
%% ============================================================

excit_p = kwave_P_time(:, 1:8:end).';
%% ============================================================
%%                   MEDIUM
%% ============================================================                  
medium.c    = 1500;
medium.rho  = 1000;
medium.beta = 0;
medium.ca   = 0;
medium.cb   = 2.0;

medium.NRL_gamma = 0.9;
medium.NRL_alpha = 0.009;

%% ============================================================
%%                   SIMULATION
%% ============================================================    

sensor_mask = zeros(mgrid.num_x, mgrid.num_y+1);
sensor_mask(:, 1:end) = 1;

tic
p_time = Backward2D(mgrid, medium, excit_p, sensor_mask, 0, 'NRL');
toc
%% ============================================================
%%                   RESULTS
%% ============================================================    
p_time = reshape(p_time, mgrid.num_t,mgrid.num_x, mgrid.num_y+1);
PAT_t(:,:) = p_time(1, :,:);

figure
subplot(1,3,3)
imagesc (mgrid.y*1e3 - min(mgrid.y*1e3), mgrid.x*1e3, (PAT_t))
axis image
colormap(jet);
caxis ([-0.05, 0.13])
colorbar 


