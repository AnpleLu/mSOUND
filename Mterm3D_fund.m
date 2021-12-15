function [M_linear, Ktemp, Ktemp2, cw] = Mterm3D_fund(mgrid, medium, omega_c)

% DESCRIPTION:
% M term calculation for 3D fundamental pressure simulation
% frequency

% USAGE:
% [M_linear, Ktemp, Ktemp2, cw] = Mterm3D_fund(mgrid, medium, omega_c)

% INPUTS:
% mgrid        input structure to define the computational domain
% medium       medium properties
% omega_c      center frequency 2*pi*fc

% OUTPUTS:
% M_linear     linear parts in M term
% Ktemp        k_z^2
% cw           frequency-dependent speed @ omegac 

%% 
if max(max(max(medium.rho))) ~= min(min(min(medium.rho)))
    % laplacian operator of density 
    rho_grad = laplacian3D(1./sqrt(medium.rho),mgrid.dx,mgrid.dy, mgrid.dz);
    rho11 = sqrt(medium.rho).*rho_grad;
    clear  rho_grad
else
    rho11 = 0;
end  
    
% transform unit from dB/cm is Np/m
alpha_Np  = abs(medium.ca.*(omega_c/2/pi/1e6).^medium.cb)*100*0.1151; 

% transform from dB/(MHz^y cm) to Np/((rad/s)^y m)
ca_Np     = 100*medium.ca.*(1e-6/(2*pi)).^medium.cb./(20*log10(exp(1))); 

% dispersion with Kramer-Kronig relation
cw        = 1./(1./medium.c + ca_Np.*tan(medium.cb*pi/2).*(abs(omega_c).^(medium.cb - 1)));   
clear ca_Np 

% diffusivity with delta = 2*cw^3*alpha/w^2
M_linear  = (medium.c0.^2./cw.^2-1).*-omega_c.^2./medium.c0^2 +...
            abs(2*alpha_Np.*(cw.^3)./(omega_c.^2))./cw.^4.*1i.*omega_c.^3 +...
            rho11;
        
clear   alpha_Np
Ktemp     = (omega_c^2./medium.c0^2) -...
            (mgrid.kx.'*ones(1, mgrid.num_y)).^2 -...
            (ones(mgrid.num_x, 1)*mgrid.ky).^2;  

kxy = squeeze((mgrid.kx.'*ones(1, mgrid.num_y)).^2 +...
             (ones(mgrid.num_x, 1)*mgrid.ky).^2);
        
Ktemp2    = (omega_c^2./cw.^2) - repmat(kxy, 1,1,mgrid.num_z+1);          
            
% Ktemp2   = (omegac^2./cw.^2) -...
%            (mgrid.kx.'*ones(1, mgrid.num_y)).^2 -...
%            (ones(mgrid.num_x, 1)*mgrid.ky).^2; 
end