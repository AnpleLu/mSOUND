function [ref_fundamental2] = BMReflection3D_fund(mgrid, medium, p_ref, ...
                             M_linear, K, K2, cw, omega_c,...
                             reflection_order, excit_p, MDMS)

% DESCRIPTION:
% Simpson-like integration for the reflection projection with correction at the
% fundamental frequency

% USAGE:
% [ref_fundamental2] = BMReflection3D_fund(mgrid, medium, p_ref, ...
%                      M_linear, K, K2, cw, omega_c,...
%                      reflection_order, excit_p, MDMS)

% INPUTS:
% mgrid             structure to define the computational domain
% medium            medium properties
% p_ref             stored reflected part
% M_linear          linear M_term
% K                 wave number along the y direction with a constant c
% K2                wave number along the y direction with a heterogeneous c
% cw                frequency-dependent speed of sound
% omega_c           transducer center frequency  
% reflection_order  the maximum order of reflection included in the simulation

% OUTPUTS:
% P_fundamental      fundamental pressure dirtribution through the domain

%%
% exponential term for the simpsion iteration
expn2 = exp(0.5*1i*K*mgrid.dz); 

% 1.05 is an empirical number, it can be higher if necessay
evanescent_fiter = 1.05;

% store all the relfected wave and its propagation
ref_fundamental = zeros(mgrid.num_x,mgrid.num_y, mgrid.num_z+1); 
ref_fundamental2 = zeros(mgrid.num_x,mgrid.num_y, mgrid.num_z+1); 

if length(cw)==1 
    cw = cw*ones(size(M_linear));
end

% c_c = medium.c;
% if length(medium.c) == 1
%     c_c = medium.c*ones(mgrid.num_x, mgrid.num_y, mgrid.num_z+1);
% end

rho_rho = medium.rho;
if length(medium.rho) == 1
    rho_rho = medium.rho*ones(mgrid.num_x, mgrid.num_y, mgrid.num_z+1);
end

h = waitbar(0,'Reflection projection, please wait...');

for iref = 1:reflection_order
    
    if mod(iref,2)==1
        p_ref2 = zeros(mgrid.num_x,mgrid.num_y, mgrid.num_z+1); 
        f = 0;
        for I = mgrid.num_z :-1:1
            
            f = f + p_ref(:,:, I+1);    
            excit_F = fftshift(fft2(f));

            Ktemp1 = (omega_c./cw(:,:,I)).^2 -...
                      evanescent_fiter*(mgrid.kx.'*ones(1, mgrid.num_y)).^2 -...
                      evanescent_fiter*(ones(mgrid.num_x, 1)*mgrid.ky).^2;
 
            % Simpson  
            % 1    
            M  = fftshift(fft2(M_linear(:,:,I+1).*f));
            F1 = excit_F.*expn2 + 0.5*mgrid.dz*M.*expn2./(2i.*K);  % P01(z+mgrid.dz/2)
            F1(isnan(F1)) = 0;
            F1(real(Ktemp1)<=0) = 0; 
            % 2
            f  = ifft2(ifftshift(F1));  % p0(z+mgrid.dz/2)
            clear F1
            M1 = fftshift(fft2(M_linear(:,:,I+1).*f));
            F2 = excit_F.*expn2 + 0.25*mgrid.dz*expn2./(2i.*K).*(M + M1./expn2);% P12(z+mgrid.dz/2)
            clear M1
            F2(isnan(F2)) = 0;
            F2(real(Ktemp1)<=0) = 0;     
            % 3   
            f  = ifft2(ifftshift(F2)); % p1(z+mgrid.dz/2)  
            M2 = fftshift(fft2(M_linear(:,:,I+1).*f));
            F3 = F2.*expn2 + 0.5*mgrid.dz*M2.*expn2./(2i.*K);  % P03(z+mgrid.dz)
            F3(isnan(F3)) = 0;
            F3(real(Ktemp1)<=0) = 0;     
            % 4     
            f  = ifft2(ifftshift(F3)); % p0(z+mgrid.dz)  
            M3 = fftshift(fft2(M_linear(:,:,I).*f));
            F4 = F2.*expn2 + 0.25*mgrid.dz.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z+mgrid.dz)
            clear M3   F3   F2
            F4(isnan(F4)) = 0;
            F4(real(Ktemp1)<=0) = 0;     
            % 5   
            f  = ifft2(ifftshift(F4)); % p0(z+mgrid.dz)  
            M4 = fftshift(fft2(M_linear(:,:,I).*f)); 
            % P25(z+mgrid.dz)
            F5 = excit_F.*exp(1i*K*mgrid.dz) +...
                 mgrid.dz/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*mgrid.dz)).*...
                 exp(1i*K*mgrid.dz)./(2i.*K); 
            clear M   M2   M4  
            F5(isnan(F5)) = 0;
            F5(real(Ktemp1)<=0) = 0; 
    
%% corrections  
            c1 = cw(:,:,I+1);
            c2 = cw(:,:,I);
            rho1 = rho_rho(:,:,I+1);
            rho2 = rho_rho(:,:,I);
               
            wx2 = omega_c.*omega_c;
            e1  = exp(1i*(K - K2(:,:,I))*mgrid.dz);
            k_dif2 = (wx2./medium.c0./medium.c0 - wx2./c2./c2)*mgrid.dz;      
            k_dif  = (wx2./medium.c0./medium.c0 - wx2./c1./c1)*mgrid.dz; 
            clear  wx2  
            correction = excit_F.*exp(1i*K2(:,:,I)*mgrid.dz).*...
                         (1-e1.*(1+k_dif./4i./K)./(1-k_dif2./4i./K));
            F5 = F5 + correction;  
            F5(isnan(F5)) = 0;
            F5(real(Ktemp1)<=0) = 0; 
         
            T_rhoc = 2.*rho2.*c2./(rho1.*c1 + rho2.*c2);  %%layered medium
            f = (ifft2(ifftshift(F5)));
            f = f.*T_rhoc;   
            
            clear F5
            
            ref_fundamental(:,:, I) = ref_fundamental(:,:,I) + f;  
%% add reflection                
            p_ref2(:,:, I+1) = (ifft2(ifftshift(excit_F))).*(T_rhoc-1); 
        end
        
    else
        p_ref = zeros(size(p_ref));
        f = 0;
        
        for I = 2:mgrid.num_z+1
            
            f = f + p_ref2(:,:,I-1); 
            excit_F = fftshift(fft2(f));           
            
            Ktemp1 = (omega_c./cw(:,:,I)).^2 -...
                      evanescent_fiter*(mgrid.kx.'*ones(1, mgrid.num_y)).^2 -...
                      evanescent_fiter*(ones(mgrid.num_x, 1)*mgrid.ky).^2;
            % Simpson  
            % 1    
            M  = fftshift(fft2(M_linear(:,:,I-1).*f));
            F1 = excit_F.*expn2 + 0.5*mgrid.dz*M.*expn2./(2i.*K);  % P01(z+mgrid.dz/2)
            F1(isnan(F1)) = 0;
            F1(real(Ktemp1)<=0) = 0; 
            % 2
            f  = ifft2(ifftshift(F1));  % p0(z+mgrid.dz/2)
            clear F1
            M1 = fftshift(fft2(M_linear(:,:,I-1).*f));
            F2 = excit_F.*expn2 + 0.25*mgrid.dz*expn2./(2i.*K).*(M + M1./expn2);% P12(z+mgrid.dz/2)
            clear M1
            F2(isnan(F2)) = 0;
            F2(real(Ktemp1)<=0) = 0;     
            % 3   
            f  = ifft2(ifftshift(F2)); % p1(z+mgrid.dz/2)  
            M2 = fftshift(fft2(M_linear(:,:,I-1).*f));
            F3 = F2.*expn2 + 0.5*mgrid.dz*M2.*expn2./(2i.*K);  % P03(z+mgrid.dz)
            F3(isnan(F3)) = 0;
            F3(real(Ktemp1)<=0) = 0;     
            % 4     
            f  = ifft2(ifftshift(F3)); % p0(z+mgrid.dz)  
            M3 = fftshift(fft2(M_linear(:,:,I).*f));
            F4 = F2.*expn2 + 0.25*mgrid.dz.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z+mgrid.dz)
            clear M3   F3   F2
            F4(isnan(F4)) = 0;
            F4(real(Ktemp1)<=0) = 0;     
            % 5   
            f  = ifft2(ifftshift(F4)); % p0(z+mgrid.dz)  
            M4 = fftshift(fft2(M_linear(:,:,I).*f)); 
            % P25(z+mgrid.dz)
            F5 = excit_F.*exp(1i*K*mgrid.dz) + ...
                 mgrid.dz/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*mgrid.dz)).*...
                 exp(1i*K*mgrid.dz)./(2i.*K); 
            clear M   M2   M4  
            F5(isnan(F5)) = 0;
            F5(real(Ktemp1)<=0) = 0; 
            
%% corrections  

            c1 = cw(:,:,I-1);
            c2 = cw(:,:,I);
            rho1 = rho_rho(:,:,I-1);
            rho2 = rho_rho(:,:,I);
               
            wx2 = omega_c.*omega_c;
            e1  = exp(1i*(K - K2(:,:,I))*mgrid.dz);
            k_dif2 = (wx2./medium.c0./medium.c0 - wx2./c2./c2)*mgrid.dz;      
            k_dif  = (wx2./medium.c0./medium.c0 - wx2./c1./c1)*mgrid.dz; 
            clear  wx2  
            correction = excit_F.*exp(1i*K2(:,:,I)*mgrid.dz).*...
                         (1-e1.*(1+k_dif./4i./K)./(1-k_dif2./4i./K));
            F5 = F5 + correction;  
            F5(isnan(F5)) = 0;
            F5(real(Ktemp1)<=0) = 0; 
         
            T_rhoc = 2.*rho2.*c2./(rho1.*c1 + rho2.*c2);  %%layered medium
            f = (ifft2(ifftshift(F5)));
            f = f.*T_rhoc;                
            clear F5

%% add reflection 
            p_ref(:,:, I-1) = (ifft2(ifftshift(excit_F))).*(T_rhoc-1); 
        end

        waitbar(iref/reflection_order)
    end
end
close (h)


if MDMS ==1
    
    backward_dz = -mgrid.dz;
    f = excit_p;
    clear excit_p
    % exponential term for the forward propagation
    expn2 = exp(0.5*1i*K*backward_dz); 

    for I = mgrid.num_z:-1:1

        rho2 = rho_rho(:,:, I);
        c1 = cw(:,:, I+1);
        c2 = cw(:,:,I);
        rho1 = rho_rho(:,:,I+1);            
        T_rhoc = 2.*rho1.*c1./(rho1.*c1 + rho2.*c2);  %%layered medium 
       
        if  (max(max(abs(T_rhoc))) ~=1 || min(min(abs(T_rhoc))) ~=1) 
            f = f - p_ref2(:,:,I);       
        end
        excit_F = fftshift(fft2(f));
    
        Ktemp1 = (omega_c./cw(:,:,I)).^2 -...
                  evanescent_fiter*(mgrid.kx.'*ones(1, mgrid.num_y)).^2 -...
                  evanescent_fiter*(ones(mgrid.num_x, 1)*mgrid.ky).^2;
        % Simpson  
        % 1    
        M  = fftshift(fft2(M_linear(:,:,I+1).*f));
        F1 = excit_F.*expn2 + 0.5*backward_dz*M.*expn2./(2i.*K);  % P01(z+backward_dz/2)
        F1(isnan(F1)) = 0;
        F1(real(Ktemp1)<=0) = 0; 
        % 2
        f  = ifft2(ifftshift(F1));  % p0(z+backward_dz/2)
        clear F1
        M1 = fftshift(fft2(M_linear(:,:,I+1).*f));
        F2 = excit_F.*expn2 + 0.25*backward_dz*expn2./(2i.*K).*(M + M1./expn2);% P12(z+backward_dz/2)
        clear M1
        F2(isnan(F2)) = 0;
        F2(real(Ktemp1)<=0) = 0;     
        % 3   
        f  = ifft2(ifftshift(F2)); % p1(z+backward_dz/2)  
        M2 = fftshift(fft2(M_linear(:,:,I+1).*f));
        F3 = F2.*expn2 + 0.5*backward_dz*M2.*expn2./(2i.*K);  % P03(z+backward_dz)
        F3(isnan(F3)) = 0;
        F3(real(Ktemp1)<=0) = 0;     
        % 4     
        f  = ifft2(ifftshift(F3)); % p0(z+backward_dz)  
        M3 = fftshift(fft2(M_linear(:,:,I).*f));
        F4 = F2.*expn2 + 0.25*backward_dz.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z+backward_dz)
        clear M3   F3   F2
        F4(isnan(F4)) = 0;
        F4(real(Ktemp1)<=0) = 0;     
        % 5   
        f  = ifft2(ifftshift(F4)); % p0(z+backward_dz)  
        M4 = fftshift(fft2(M_linear(:,:,I).*f)); 
        % P25(z+backward_dz)
        F5 = excit_F.*exp(1i*K*backward_dz) + ...
             backward_dz/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*backward_dz)).*...
             exp(1i*K*backward_dz)./(2i.*K); 
        clear M   M2   M4  
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
    
%% corrections  
               
        wx2 = omega_c.*omega_c;
        e1  = exp(1i*(K - K2(:,:,I))*backward_dz);
        k_dif2 = (wx2./medium.c0./medium.c0 - wx2./c2./c2)*backward_dz;      
        k_dif  = (wx2./medium.c0./medium.c0 - wx2./c1./c1)*backward_dz; 
        clear  wx2  
        correction = excit_F.*exp(1i*K2(:,:,I)*backward_dz).*...
                    (1-e1.*(1+k_dif./4i./K)./(1-k_dif2./4i./K));
        F5 = F5 + correction;  
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
        
        f = (ifft2(ifftshift(F5)));
        f = f./T_rhoc;                
                
        F5 = fftshift(fft2(f));
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 

%         excit_F = F5;
        f = ifft2(ifftshift(F5));
        clear F5  
        
        % recover pressure from the normalized wave field
        ref_fundamental2(:,:,I) = ref_fundamental(:,:,I) + f;
    end
else 
    ref_fundamental2 = ref_fundamental;        
end

end
