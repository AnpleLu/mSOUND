function [M_nonlinear, M_linear, cw] = Mterm3D(mgrid, medium, medium_cond)
% DESCRIPTION:
% calculate the M term at RHS of solution equation

% USAGE:
% [M_nonlinear, M_linear, cw] = Mterm3D(mgrid, medium, medium_cond)

% INPUTS:
% mgrid    
% medium       medium properties
% medium_cond  returned from function medium_case 

% OUTPUTS:
% M_lienar     linear terms in the M
% M_nonlinear  nonlinear term in the M
% cw           speed of sound when dispersion is considered
%%
if medium_cond == 1
    
    M_nonlinear = 0;
    M_linear   = single(zeros(mgrid.num_t, 1));
    medium.ca0 = medium.ca(1);
    medium.cb0 = medium.cb(1);   
    
    % transform from dB/(MHz^y cm) to Np/((rad/s)^y m)
    ca_Np = 100*medium.ca0.*(1e-6/(2*pi)).^medium.cb0./(20*log10(exp(1)));
    
    % calculate frequency-dependent speed with Kramers-Kronig relation
    cw = 1./(1./medium.c0 + ca_Np.*tan(medium.cb0*pi/2)...
         .*(abs(mgrid.w.').^(medium.cb0 - 1) - 0));      
    
    for it = 1:mgrid.num_t  
        % unit is dB/cm and then transformed to Np/m
        alpha_Np = abs(medium.ca0.*((mgrid.w(it))/2/pi/1e6).^medium.cb0)*100*0.1151; 
        M_linear(it,1) = abs(2*alpha_Np.*(medium.c0.^3)./((mgrid.w(it)).^2))./...
                         medium.c0.^4.*1i.*(mgrid.w(it)).^3 ;
        clear  alpha_Np               
    end

 %%
elseif medium_cond == 2 % 'homogeneous medium'

    medium.rho0 = medium.rho(1);    % Density of medium [kg/m^3]
    medium.beta0 = medium.beta(1);
    medium.ca0 = medium.ca(1);
    medium.cb0 = medium.cb(1);
    
    M_nonlinear = zeros(mgrid.num_t, 1);
    M_linear    = zeros(mgrid.num_t, 1);
    cw          = zeros(mgrid.num_t, 1);
    
    % transform from dB/(MHz^y cm) to Np/((rad/s)^y m)
    ca_Np = 100*medium.ca0.*(1e-6/(2*pi)).^medium.cb0./(20*log10(exp(1)));
     
    for it = 1:mgrid.num_t
    
        wt = mgrid.w(it);
        % unit is dB/cm and then transformed to Np/m
    
        alpha_Np = abs(medium.ca0.*(wt/2/pi/1e6).^medium.cb0)*100*0.1151; 
    
        % calculate frequency-dependent speed with Kramers-Kronig relation
        cwt      = 1./(1./medium.c0 + ca_Np.*tan(medium.cb0*pi/2)...
                  .*(abs(wt).^(medium.cb0 - 1) - 0)); 
           
        M_nonlinear(it, 1) = medium.beta0./(medium.rho0)./cwt.^4.*wt.^2;  
    
        % diffusivity with delta = 2*cw^3*alpha/w^2       
        M_linear(it,1) = abs(2*alpha_Np.*(cwt.^3)./(wt.^2))./cwt.^4.*1i.*wt.^3 ;
        clear wt  alpha_Np               
                           
        cw(it,1) = cwt;
     end
     clear  ca_Np   wt   alpha_Np  cwt      
    
elseif medium_cond == 3 % 'only nonlinear coefficient variation'

    medium.rho0 = medium.rho(1);    % Density of medium [kg/m^3]
    medium.beta0 = medium.beta(1);
    medium.ca0 = medium.ca(1);
    medium.cb0 = medium.cb(1);   
    
    M_nonlinear = zeros(mgrid.num_t, mgrid.num_x, mgrid.num_y, mgrid.num_z+1);
    M_linear    = zeros(mgrid.num_t, 1);
    cw          = zeros(mgrid.num_t, 1);

    % transform from dB/(MHz^y cm) to Np/((rad/s)^y m)
    ca_Np = 100*medium.ca0.*(1e-6/(2*pi)).^medium.cb0./(20*log10(exp(1)));  

    for it = 1:mgrid.num_t
    
        wt = mgrid.w(it);

        % calculate frequency-dependent speed with Kramers-Kronig relation
        cwt = 1./(1./medium.c0 + ca_Np.*tan(medium.cb0*pi/2)...
                  .*(abs(wt).^(medium.cb0 - 1) - 0)); 
           
        % unit is dB/cm and then transformed to Np/m
        alpha_Np = abs(medium.ca0.*(wt/2/pi/1e6).^medium.cb0)*100*0.1151; 
        
        % diffusivity with delta = 2*cw^3*alpha/w^2       
        M_linear(it,:) = abs(2*alpha_Np.*(cwt.^3)./(wt.^2))./cwt.^4.*1i.*wt.^3 ;

        M_nonlinear(it, :, :,:) = medium.beta./(medium.rho0)./cwt.^4.*wt.^2;                     
        cw(it,:) = cwt;
        clear wt  alpha_Np  
     end
     clear  ca_Np   wt   alpha_Np  cwt    
    
elseif medium_cond == 4  % 'linear media'

    if max(max(max(medium.rho))) ~= min(min(min(medium.rho)))
        % laplacian operator of density 
        rho_grad = laplacian3D(1./sqrt(medium.rho),mgrid.dx,mgrid.dy, mgrid.dz);
        rho11 = sqrt(medium.rho).*rho_grad;
        clear  rho_grad
    else
        rho11 = 0;
    end      

    M_nonlinear = 0;
    M_linear = zeros(mgrid.num_t, mgrid.num_x, mgrid.num_y, mgrid.num_z+1);
    cw = zeros(mgrid.num_t, mgrid.num_x, mgrid.num_y, mgrid.num_z+1);

    % transform from dB/(MHz^y cm) to Np/((rad/s)^y m)
    ca_Np = 100*medium.ca.*(1e-6/(2*pi)).^medium.cb./(20*log10(exp(1)));  
    
    for it = 1:mgrid.num_t
    
        wt = mgrid.w(it)*ones(mgrid.num_x, mgrid.num_y, mgrid.num_z+1);
        % unit is dB/cm and then transformed to Np/m
    
        alpha_Np = abs(medium.ca.*(wt/2/pi/1e6).^medium.cb)*100*0.1151; 
    
        % calculate frequency-dependent speed with Kramers-Kronig relation
        cwt      = 1./(1./medium.c + ca_Np.*tan(medium.cb*pi/2)...
               .*(abs(wt).^(medium.cb - 1) - 0)); 
           
        % diffusivity with delta = 2*cw^3*alpha/w^2       
        M_linear(it,:,:,:) = (medium.c0.^2./cwt.^2-1).*-wt.^2./medium.c0.^2 +...
                             abs(2*alpha_Np.*(cwt.^3)./(wt.^2))./cwt.^4.*1i.*wt.^3 + ... 
                             rho11;
        cw(it,:,:,:) = cwt;
        clear wt  alpha_Np cwt   
     end    
    
 %%                   
elseif medium_cond == 5 % 'nonlinear media'

    if max(max(max(medium.rho))) ~= min(min(min(medium.rho)))
        % laplacian operator of density 
        rho_grad = laplacian3D(1./sqrt(medium.rho),mgrid.dx, mgrid.dy, mgrid.dz);
        rho11 = sqrt(medium.rho).*rho_grad;
        clear  rho_grad
    else
        rho11 = 0;
    end      

    M_nonlinear = single(zeros(mgrid.num_t, mgrid.num_x, mgrid.num_y, mgrid.num_z+1));
    cw = zeros(mgrid.num_t, mgrid.num_x, mgrid.num_y, mgrid.num_z+1);
    M_linear = zeros(mgrid.num_t, mgrid.num_x, mgrid.num_y, mgrid.num_z+1);

    % transform from dB/(MHz^y cm) to Np/((rad/s)^y m)
    ca_Np = 100*medium.ca.*(1e-6/(2*pi)).^medium.cb./(20*log10(exp(1)));  

    for it = 1:mgrid.num_t
    
        wt = mgrid.w(it)*ones(mgrid.num_x, mgrid.num_y, mgrid.num_z+1);
        % unit is dB/cm and then transformed to Np/m
    
        alpha_Np = abs(medium.ca.*(wt/2/pi/1e6).^medium.cb)*100*0.1151; 
    
        % calculate frequency-dependent speed with Kramers-Kronig relation
        cwt      = 1./(1./medium.c + ca_Np.*tan(medium.cb*pi/2)...
                  .*(abs(wt).^(medium.cb - 1) - 0)); 
           
        M_nonlinear(it, :, :,:) = medium.beta./sqrt(medium.rho)./cwt.^4.*wt.^2;  
    
        % diffusivity with delta = 2*cw^3*alpha/w^2       
        M_linear(it,:,:,:) = (medium.c0.^2./cwt.^2-1).*-wt.^2./medium.c0.^2 +...
                             abs(2*alpha_Np.*(cwt.^3)./(wt.^2))./cwt.^4.*1i.*wt.^3 + ... 
                             rho11;
        clear wt  alpha_Np               
                           
        cw(it,:,:,:) = cwt;
     end
     clear  ca_Np   wt   alpha_Np  cwt 

end


             
end




