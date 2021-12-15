clear
% Simulation of a strongly 2D heterogeneous medium 
% using the transient mixed domain method
%% ============================================================
%%                   PARAMETERS
%% ============================================================

fc = 0.7e6;	         % fundamental frequency [Hz]
p0 = 1e6;	         % initial pressure amplitude [Pa]
medium.c0 = 1500;    % speed of Sound [m/s]
medium.rho0 = 1000;  % density of medium [kg/m^3]

omegac = fc*2*pi;      % fundamental angular frequency
lambda = medium.c0/fc; % wavelength [m]

dx = lambda/4;     % step size in the x direction [m]
dy = lambda/16;    % step size in the y direction [m]
dt = 1/fc/12;      % timestep [s]

x_length = lambda*78; % computational domain size in the x direction [m]
y_length = lambda*40; % computational domain size in the y direction [m]
t_length = 4/fc*15;   % temporal domain size [s]

%% ============================================================
%%                   COMPUTATIONAL DOMAIN
%% ============================================================

mgrid = set_grid(dt, t_length, dx, x_length, dy, y_length);

%% ============================================================
%%                   TRANSDUCER
%% ============================================================

TR_focus  = 30*lambda;    % transducer focal length [m]
TR_radius = 10*lambda;    % transducer radius [m]
focus  = round(TR_focus/mgrid.dy);       % [grid points]
num_el = round(2*TR_radius/mgrid.dx);    % [grid points]

el_x = -floor(num_el/2):floor(num_el/2); % [grid points]
delay = sqrt((el_x*mgrid.dx).^2 + (focus*mgrid.dy)^2)/medium.c0;  % [s]
delay = delay - min(delay);

%% ============================================================
%%                   EXCITATION
%% ============================================================

ts = [-4/fc:dt:4/fc].';
excit_ps = zeros(length(ts),num_el);
for cell = 1:num_el
    excit_ps(:,cell) = p0*sin(2*pi*fc*(ts+delay(cell))).*exp(-(ts+delay(cell)).^2*fc^2/2);
end

%% ============================================================
%%                   MEDIA
%% ============================================================

load 2D_strongly_heterogeneous_model.mat
medium.c    = c;
medium.rho  = rho;
medium.beta = beta;
medium.ca   = ca;       % power law relation constant     
medium.cb   = cb;

figure
imagesc(mgrid.y*1e3 - min(mgrid.y)*1e3, mgrid.x*1e3,medium.c)
colormap (jet)
colorbar
axis equal
axis image
hold on
plot((TR_focus)*1000,0,'r.','MarkerSize',20)
hold off
graph1 = line([0,0],[-TR_radius*1000,TR_radius*1000], 'Color','red');
set(graph1,'LineWidth',4);
axis equal

%% ============================================================
%%                   SIMULATION
%% ============================================================
sensor_mask = zeros(mgrid.num_x, mgrid.num_y+1);
sensor_mask(81:end-80, 2:end) = 1;

% the maximum order of reflection included in the simulatio
reflection_order = 4;

% 2D forward simulation  
p_time = Forward2D(mgrid, medium, excit_ps, sensor_mask, reflection_order, 'correction');

%% ============================================================
%%                   RESULTS
%% ============================================================
p_time = reshape(p_time, mgrid.num_t,152, mgrid.num_y);
P_freq = zeros(size(p_time));
for inumz = 1:mgrid.num_y
    for inumx = 1:152
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
imagesc(mgrid.y*1e3, -mgrid.x(81:end-80)*1e3, P_h1_ASA)
xlabel ('y (mm)')
ylabel ('x (mm)')
colormap(jet)

axis image
subplot(1,2,2)
imagesc(mgrid.y*1e3, -mgrid.x(81:end-80)*1e3, P_h2_ASA)
xlabel ('y (mm)')
ylabel ('x (mm)')
colormap(jet)
axis image
  
focus = round(TR_focus/mgrid.dy);
P_time(:) = p_time(:, round(152/2), focus);
figure
subplot(1,2,1)
plot(mgrid.t*1e6, P_time/1e6, 'k', 'LineWidth',2)
xlabel ('Time (\mus)')
ylabel ('Pressure (MPa)')

subplot(1,2,2)
plot(mgrid.w/2/pi/1e6, 20*log10(abs(fftshift(fftn(P_time)))), 'k', 'LineWidth',2)
xlabel ('Frequency (MHz)')
ylabel ('Pressure (dB)')




