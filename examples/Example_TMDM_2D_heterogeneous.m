clear
% Simulation of a 2D heterogeneous medium 
% using the transient mixed domain method
%% ======================================================================== 
%%                        PARAMETERS 
%% ========================================================================
fc = 1.0e6;	           % fundamental frequency	[Hz]
p0 = 1.0e6;	           % initial pressure amplitude  [Pa]
medium.c0 = 1500;      % speed of Sound [m/s]
medium.rho0 = 1000;    % density of medium [kg/m^3]

omegac = fc*2*pi;      % fundamental angular frequency
lambda = medium.c0/fc; % wave length [m]

dx = lambda/8;         % step size in the x direction [m]
dy = lambda/8;         % step size in the y direction [m]
dt = 1/fc/16;          % timestep  [s]

x_length = lambda*50;  % computational domain size in the x direction [m]
y_length = lambda*22;  % computational domain size in the y direction [m]
t_length = 4/fc*15.0;  % temporal domain size [s]

%% ======================================================================== 
%%                        COMPUTATIONAL DOMAIN
%% ========================================================================

mgrid = set_grid(dt, t_length, dx, x_length, dy, y_length);  %% wavevector and x, y

%% ======================================================================== 
%%                        TRANSDUCER
%% ========================================================================

TR_focus  = 14*lambda;     % transducer focal length [m]
TR_radius = 7.5*lambda;    % Transducer radius [m]

% calculate the phase delay
delay = sqrt((mgrid.x).^2 + (TR_focus)^2)/medium.c0;  % [s]
delay = delay - min(delay);

%% ======================================================================== 
%%                        EXCITATION
%% ========================================================================

% define the pulse length
ts = [-4/fc:dt:4/fc].';
delay = repmat(delay, length(ts),1);
ts = repmat(ts, 1,mgrid.num_x);

% define the excitation pulse
source_p = p0*sin(2*pi*fc*(ts+delay)).*exp(-(ts+delay).^2*fc^2/2);
source_p(:, abs(mgrid.x)>TR_radius) = 0;
%% ======================================================================== 
%%                        MEDIUM
%% ========================================================================
load tissue_model.mat
medium.c    = c;
medium.rho  = rho;
medium.beta = beta;
medium.ca   = ca;       % power law relation constant     
medium.cb   = cb;
medium.NRL_gamma = 0.5;
medium.NRL_alpha = 0.05;

% visualize the 2D model
figure
imagesc(mgrid.y*1e3, mgrid.x*1e3,beta)
colormap(jet)
axis image
hold on
plot((TR_focus+min(mgrid.y))*1000,0,'r.','MarkerSize',20)
graph1 = line([min(mgrid.y)*1e3,min(mgrid.y)*1e3],[-TR_radius*1000,TR_radius*1000], 'Color','red');
set(graph1,'LineWidth',4);
axis equal
%% ======================================================================== 
%%                        SIMULATION
%% ========================================================================
% define where the pressure will be recorded
sensor_mask = zeros(mgrid.num_x, mgrid.num_y+1);
sensor_mask(100:300, 2:end) = 1;

% the maximum order of reflection included in the simulation
reflection_order = 2; 

% 2D forward simulation
p_time = Forward2D(mgrid, medium, source_p, sensor_mask, reflection_order, 'NRL');

%% ======================================================================== 
%%                        RESULTS
%% ========================================================================

p_time = reshape(p_time, mgrid.num_t, 201, mgrid.num_y);

P_freq = zeros(size(p_time));
for inumz = 1:mgrid.num_y
    for inumx = 1:201
        P_freq(:,inumx, inumz) = abs(fftshift(fft(p_time(:, inumx, inumz))));
    end
end
f1 = find(abs(-mgrid.w/2/pi-fc)   == min(abs(-mgrid.w/2/pi-fc)));
f2 = find(abs(-mgrid.w/2/pi-2*fc) == min(abs(-mgrid.w/2/pi-2*fc)));
P_h1_ASA(:,:) = P_freq(f1,:,:);
P_h2_ASA(:,:) = P_freq(f2,:,:);

% visualize the results
figure
subplot(1,2,1)
imagesc(mgrid.y*1e3, -mgrid.x(100:300)*1e3, P_h1_ASA)
xlabel ('y (mm)')
ylabel ('x (mm)')
colormap(jet)

axis image
subplot(1,2,2)
imagesc(mgrid.y*1e3, -mgrid.x(100:300)*1e3, P_h2_ASA)
xlabel ('y (mm)')
ylabel ('x (mm)')
colormap(jet)
axis image

focus = round(TR_focus/mgrid.dy);
P_time(:) = p_time(:, round(201/2), focus);
figure
subplot(1,2,1)
plot(mgrid.t*1e6, P_time/1e6, 'k', 'LineWidth',2)
xlabel ('Time (\mus)')
ylabel ('Pressure (MPa)')

subplot(1,2,2)
plot(mgrid.w/2/pi/1e6, 20*log10(abs(fftshift(fftn(P_time)))), 'k', 'LineWidth',2)
xlabel ('Frequency (MHz)')
ylabel ('Pressure (dB)')

