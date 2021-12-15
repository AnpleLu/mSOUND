clear
% Simulation of a 2D homogeneous medium 
% using the transient mixed domain method
%% ======================================================================== 
%%                        PARAMETERS 
%% ========================================================================
fc = 1.0e6;	           % fundamental frequency	[Hz]
p0 = 1;	               % initial pressure amplitude [Pa]
medium.c0 = 1500;      % speed of Sound [m/s]
medium.rho0 = 1000;    % density of medium [kg/m^3]

omegac = fc*2*pi;      % fundamental angular frequency
lambda = medium.c0/fc; % wavelength [m]

dx = lambda/4;         % step size in the x direction [m] 
dy = lambda/4;         % step size in the y direction [m] 
dt = 1/fc/16;          % timestep  [s]

x_length = lambda*30+dx;  % computational domain size in the x direction [m] 
y_length = lambda*16;     % computational domain size in the y direction [m] 
t_length = 4/fc*15.0;     % temporal domain size [s] 

%% ======================================================================== 
%%                        COMPUTATIONAL DOMAIN
%% ========================================================================

mgrid = set_grid(dt, t_length, dx, x_length, dy, y_length); 

%% ======================================================================== 
%%                        TRANSDUCER
%% ========================================================================

TR_focus  = 12*lambda;  % transducer focal length [m]
TR_radius = 6*lambda;   % Transducer radius [m]
focus  = round(TR_focus/mgrid.dy);   % [grid points]

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

medium.c    = 1500;   % speed of sound [m/s] 
medium.rho  = 1000;   % density [kg/m^3] 
medium.beta = 0;      % nonlinearity coefficient
medium.ca   = 0.0;    % power law relation constant     
medium.cb   = 2.0;    % power law exponent

% non-reflecting layer
medium.NRL_gamma = 0.5;
medium.NRL_alpha = 0.05;
%% ======================================================================== 
%%                        SIMULATION
%% ========================================================================
% define where the pressure will be recorded
sensor_mask = zeros(mgrid.num_x, mgrid.num_y+1);
sensor_mask(30:90, 2:end) = 1;

% 2D forward simulation
p_time = Forward2D(mgrid, medium, source_p, sensor_mask, 0, 'NRL');

%% ======================================================================== 
%%                        RESULTS
%% ========================================================================

p_time = reshape(p_time, mgrid.num_t,61, mgrid.num_y);

P_freq = zeros(size(p_time));
for inumz = 1:mgrid.num_y
    for inumx = 1:61
        P_freq(:,inumx, inumz) = abs(fftshift(fft(p_time(:, inumx, inumz))));
    end
end
f1 = find(abs(-mgrid.w/2/pi-fc)   == min(abs(-mgrid.w/2/pi-fc)));
P_h1_ASA(:,:) = P_freq(f1,:,:);

% visualize the results
figure
p_focus = squeeze(p_time(:, round(61/2), end));
subplot(1,2,1)
plot(mgrid.t*1e6, p_focus, 'k', 'LineWidth',2)
xlabel ('Time (\mus)')
ylabel ('Pressure (Pa)')

subplot(1,2,2)
imagesc(mgrid.y*1e3, mgrid.x(30:90)*1e3, P_h1_ASA)
xlabel ('y (mm)')
ylabel ('x (mm)')
colormap(jet)
axis image
