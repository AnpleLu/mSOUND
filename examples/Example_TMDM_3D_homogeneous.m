clear
% Simulation of a 3D homogeneous medium 
% using the transient mixed domain method
%% ======================================================================== 
%%                        PARAMETERS 
%% ========================================================================
fc = 1.0e6;	           % fundamental frequency	[Hz]
p0 = 1.0e6;	           % initial pressure amplitude [Pa]
medium.c0 = 1500;      % speed of Sound [m/s]

omegac = fc*2*pi;      % fundamental angular frequency
lambda = medium.c0/fc; % wavelength [m]

dx = lambda/8;  % step size in the x direction [m]
dy = lambda/8;  % step size in the y direction [m]
dz = lambda/8;  % step size in the z direction [m]
dt = 1/fc/16;   % timestep [s]

x_length = lambda*20; % computational domain size in the x direction [m]
y_length = lambda*20; % computational domain size in the y direction [m]
z_length = lambda*10; % computational domain size in the z direction [m]
t_length = 4/fc*5.0;  % temporal domain size [s] 

%% ======================================================================== 
%%                        COMPUTATIONAL DOMAIN
%% ========================================================================

mgrid = set_grid(dt, t_length, dx, x_length, dy, y_length, dz, z_length);  

%% ======================================================================== 
%%                        TRANSDUCER
%% ========================================================================
TR_focus  = 8*lambda;               % transducer focal length [m]
TR_radius = 5.0*lambda;             % transducer radius [m]
focus  = round(TR_focus/mgrid.dy);  % [grid points]

% calculate the phase delay
X = ones(mgrid.num_y,1)*mgrid.x; 
Y = (ones(mgrid.num_x,1)*mgrid.y).'; 
delay = sqrt(X.^2 + Y.^2 + TR_focus^2)/medium.c0; % calculate the delay for each element [m] 
delay = delay - min(min(delay)); 

%% ======================================================================== 
%%                        EXCITATION
%% ========================================================================
ts = [-4/fc:dt:4/fc].';
source_p = zeros(length(ts), mgrid.num_x, mgrid.num_y);
for cel1 = 1:mgrid.num_x
    for cel2 = 1:mgrid.num_y
        source_p(:, cel1, cel2) = ...
        p0*sin(2*pi*fc*(ts + delay(cel1, cel2))).*...
        exp(-(ts+ delay(cel1, cel2)).^2*fc^2/2);
    end
end

% generate a circular transducer
RHO = sqrt(X.^2+Y.^2);

source_p(:, RHO>TR_radius) = 0;
clear excit_ptemp X Y
%% ======================================================================== 
%%                        MEDIUM
%% ========================================================================                 

medium.c    = 1500;
medium.rho  = 1000;
medium.beta = 3.6;
medium.ca   = 0.005; % power law relation constant     
medium.cb   = 2.0;

%% ======================================================================== 
%%                        SIMULATION
%% ========================================================================
% define where the pressure will be recorded
sensor_mask = zeros(mgrid.num_x, mgrid.num_y, mgrid.num_z+1);
sensor_mask(40:120, 40:120, 2:end) = 1;

% 3D forward simulation
p_total = Forward3D(mgrid, medium, source_p, sensor_mask, 0);

p_total = reshape(p_total, mgrid.num_t, 81, 81, mgrid.num_z);

%% ======================================================================== 
%%                        RESULTS
%% ========================================================================
P_freq = zeros(size(p_total));
for inumz = 1:mgrid.num_z
    for inumy = 1:81
        for inumx = 1:81
            P_freq(:,inumx, inumy, inumz) = abs(fftshift(fft(p_total(:, inumx, inumy, inumz))));
        end
    end
end
f1 = find(abs(-mgrid.w/2/pi-fc)   == min(abs(-mgrid.w/2/pi-fc)));
f2 = find(abs(-mgrid.w/2/pi-2*fc) == min(abs(-mgrid.w/2/pi-2*fc)));
P_h1_ASAxz(:,:) = P_freq(f1,round(81/2),:,:);
P_h1_ASAxy(:,:) = P_freq(f1,:,:,focus);
P_h2_ASAxz(:,:) = P_freq(f2,round(81/2),:,:);
P_h2_ASAxy(:,:) = P_freq(f2,:,:,focus);

% visualize the results
figure
subplot(2,2,1)
imagesc(mgrid.z*1e3, mgrid.x(40:120)*1e3, P_h1_ASAxz)
xlabel ('y (mm)')
ylabel ('x (mm)')
colormap(jet)
axis image
subplot(2,2,2)
imagesc(mgrid.y(40:120)*1e3, mgrid.x(40:120)*1e3,P_h1_ASAxy)
xlabel ('y (mm)')
ylabel ('x (mm)')
colormap(jet)
axis image
subplot(2,2,3)
imagesc(mgrid.z*1e3, mgrid.x(40:120)*1e3,P_h2_ASAxz)
xlabel ('y (mm)')
ylabel ('x (mm)')
colormap(jet)
axis image
subplot(2,2,4)
imagesc(mgrid.y(40:120)*1e3, mgrid.x(40:120)*1e3,P_h2_ASAxy)
xlabel ('y (mm)')
ylabel ('x (mm)')
colormap(jet)
axis image  

P_time(:) = p_total(:,round(81/2),round(81/2),end);
figure
subplot(1,2,1)
plot(mgrid.t*1e6, P_time/1e6, 'k', 'LineWidth',3)
P_dB = 20*log10(abs(fftshift(fftn(P_time))));
subplot(1,2,2)
plot(mgrid.w/2/pi/1e6, P_dB, 'k', 'LineWidth',3)
