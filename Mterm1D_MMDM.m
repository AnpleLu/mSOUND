function [M_nonlinear, M_linear, cw] = Mterm1D_MMDM(mgrid, medium, medium_cond)
% DESCRIPTION:
% calculate the M term at RHS of solution equation

% USAGE:
%  [M_nonlinear, M_linear, cw] = Mterm1D_MMDM(mgrid, medium, medium_cond)

% INPUTS:
% mgrid    
% medium       medium properties
% medium_cond  returned from function medium_case 

% OUTPUTS:
% M_lienar     linear terms in the M
% M_nonlinear  nonlinear term in the M
% cw           speed of sound when dispersion is considered

%%

if medium_cond == 4  % 'linear media'
    
    M_nonlinear = 0;
    cw          = zeros(mgrid.num_t, mgrid.num_x+1);
    M_linear    = zeros(mgrid.num_t, mgrid.num_x+1);
    % transform from dB/(MHz^y cm) to Np/((rad/s)^y m)
    ca_Np = 100*medium.ca.*(1e-6/(2*pi)).^medium.cb./(20*log10(exp(1)));  

    for it = 1:mgrid.num_t
        
        wt = mgrid.w(it)*ones(1,mgrid.num_x+1);
    
        % unit is dB/cm and then transformed to Np/m
        alpha_Np = abs(medium.ca.*(wt/2/pi/1e6).^medium.cb)*100*0.1151;
    
        % calculate frequency-dependent speed with Kramers-Kronig relation
        cwt = 1./(1./medium.c + ca_Np.*tan(medium.cb*pi/2)...
              .*(abs(wt).^(medium.cb - 1) - 0)); 
           
        % calculate the linear parts in M term 
        % diffusivity with delta = 2*cw^3*alpha/w^2
        M_linear(it,:) = (medium.c0.^2./cwt.^2-1).*-wt.^2./medium.c0.^2 +...
                            abs(2*alpha_Np.*(cwt.^3)./(wt.^2))./cwt.^4.*1i.*wt.^3;
        
        clear  alpha_Np   wt
        cw(it,:)  = cwt;
        clear cwt  
    end   
   
elseif medium_cond == 5 % 'nonlinear media'
    
    M_nonlinear = zeros(mgrid.num_t, mgrid.num_x+1);
    cw          = zeros(mgrid.num_t, mgrid.num_x+1);
    M_linear    = zeros(mgrid.num_t, mgrid.num_x+1);
    
    % transform from dB/(MHz^y cm) to Np/((rad/s)^y m)
    ca_Np = 100*medium.ca.*(1e-6/(2*pi)).^medium.cb./(20*log10(exp(1)));  

    for it = 1:mgrid.num_t
        
        wt = mgrid.w(it)*ones(1, mgrid.num_x+1);
    
        % unit is dB/cm and then transformed to Np/m
        alpha_Np = abs(medium.ca.*(wt/2/pi/1e6).^medium.cb)*100*0.1151;
    
        % calculate frequency-dependent speed with Kramers-Kronig relation
        cwt = 1./(1./medium.c + ca_Np.*tan(medium.cb*pi/2)...
              .*(abs(wt).^(medium.cb - 1) - 0)); 
          
        M_nonlinear(it, :) = medium.beta./(medium.rho)./cwt.^4.*wt.^2;
        
        % calculate the linear parts in M term 
        % diffusivity with delta = 2*cw^3*alpha/w^2
        M_linear(it,:) = (medium.c0.^2./cwt.^2-1).*-wt.^2./medium.c0.^2 +...
                           abs(2*alpha_Np.*(cwt.^3)./(wt.^2))./cwt.^4.*1i.*wt.^3;
   
        clear  alpha_Np   wt
                                
        cw(it,:)  = cwt;
        clear cwt  
    end
      
end

end
