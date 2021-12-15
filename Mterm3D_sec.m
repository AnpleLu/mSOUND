function [M_linear, Ktemp2, cw2, cw] = Mterm3D_sec(mgrid, medium, omega_c)

% DESCRIPTION:
% M term calculation for 3D second-harmonic pressure simulation
% frequency

% USAGE:
% [M_linear, Ktemp2, cw2, cw] = Mterm3D_sec(mgrid, medium, omega_c)

% INPUTS:
% mgrid        input structure to define the computational domain
% medium       medium properties
% omega_c      center frequency 2*pi*fc

% OUTPUTS:
% M_linear     linear parts in M term
% Ktemp2       k_z^2  @ 2*omegac
% cw2          frequency-dependent speed @ 2*omegac
% cw           frequency-dependent speed @ omegac
%%
omeagc2    = 2*omega_c;
Ktemp2     = ((omeagc2).^2./(medium.c0).^2)*ones(mgrid.num_x,mgrid.num_y) -...
             (mgrid.kx.'*ones(1, mgrid.num_y)).^2 - ...
             (ones(mgrid.num_x, 1)*mgrid.ky).^2;

% transform unit from dB/cm is Np/m
alpha_Np   = abs(medium.ca.*(omeagc2/2/pi/1e6).^(medium.cb))*100*0.1151;      

% transform from dB/(MHz^y cm) to Np/((rad/s)^y m)
ca_Np      = 100*(medium.ca).*(1e-6/(2*pi)).^(medium.cb)./(20*log10(exp(1))); 

% dispersion with Kramer-Kronig relation @ 2*omegac
cw2 = 1./(1./(medium.c) + ca_Np.*tan(medium.cb*pi/2)...
      .*(abs(omeagc2).^(medium.cb - 1) )); 
  
% dispersion with Kramer-Kronig relation @ omegac  
cw  = 1./(1./medium.c + ca_Np.*tan(medium.cb*pi/2)...
      .*(abs(omega_c).^(medium.cb - 1)));   
 
clear  ca_Np 
% linear term in M
M_linear   = (medium.c0.^2./cw2.^2-1).*-omeagc2.^2./medium.c0^2 +...
             1i*abs(2*alpha_Np.*(cw2.^3)./(omeagc2.^2))./cw2.^4*omeagc2.^3 ;
clear alpha_Np
end