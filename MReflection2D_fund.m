function [ref_fundamental] = MReflection2D_fund(mgrid, medium, p_ref, ...
                             M_linear, K, K2, cw, omegac, reflection_order)

% DESCRIPTION:
% Simpson integration for the reflection projection with correction at the
% fundamental frequency

% USAGE:
% [ref_fundamental] = MReflection2D_fund(mgrid, medium, p_ref, ...
%                     M_linear, K, K2, cw, omegac, reflection_order)

% INPUTS:
% mgrid              structure to define the computational domain
% medium             medium properties
% p_ref              stored reflected part
% M_linear           linear M_term
% K                  wave number along the y direction with a constant c
% K2                 wave number along the y direction with a heterogeneous c
% cw                 frequency-dependent speed of sound
% omegac             transducer center frequency  
% reflection_order   the maximum order of reflection included in the simulation


% OUTPUTS:
% P_fundamental      fundamental pressure dirtribution through the domain

%%
% exponential term for the simpsion iteration
expn2 = exp(0.5*1i*K*mgrid.dy); 

% store all the relfected wave and its propagation
ref_fundamental = zeros(mgrid.num_x,mgrid.num_y+1); 

% 1.05 is an empirical number, it can be higher if necessay
evanescent_fiter = 1.05;

c_c = medium.c;
if length(medium.c) == 1
    c_c = medium.c*ones(mgrid.num_x, mgrid.num_y+1);
end

rho_rho = medium.rho;
if length(medium.rho) == 1
    rho_rho = medium.rho*ones(mgrid.num_x, mgrid.num_y+1);
end

h = waitbar(0,'Reflection projection, please wait...');

for iref = 1:reflection_order
    
    if mod(iref,2)==1
        p_ref2 = zeros(mgrid.num_x,mgrid.num_y+1); 
        f = 0;
        for I = mgrid.num_y:-1:1 % the layer which the main waveform arrives 
            
            c1 = c_c(:,I+1);
            c2 = c_c(:,I);
            rho1 = rho_rho(:,I+1);
            rho2 = rho_rho(:,I);
    
            Ktemp = omegac.^2./cw(:,I).^2 - evanescent_fiter*(mgrid.kx.').^2;    
                
            f = f + p_ref(:,I+1);      
            excit_F = fftshift(fft(f));

            %Simpson   
            % 1    
            M  = fftshift(fft(M_linear(:,I+1).*f));           % M(f(z))
            F1 = excit_F.*expn2 + 0.5*mgrid.dy*M.*expn2./(2i.*K);  % P01(z-mgrid.dy/2)
            F1(isnan(F1)) = 0;
            F1(real(Ktemp)<=0) = 0; 

            % 2
            f  = (ifft(ifftshift(F1)));  % p0(z-mgrid.dy/2)
            M1 = fftshift(fft(M_linear(:,I+1).*f));  % M(f(z))
            F2 = excit_F.*expn2 + 0.25*mgrid.dy*expn2./(2i.*K).*(M + M1./expn2);% P12(z-mgrid.dy/2)
            F2(isnan(F2)) = 0;
            F2(real(Ktemp)<=0) = 0; 
    
            % 3   
            f  = (ifft(ifftshift(F2))); % p1(z-mgrid.dy/2)  
            M2 = fftshift(fft(M_linear(:,I+1).*f));  % M(f1(z-mgrid.dy/2))
            F3 = F2.*expn2 + 0.5*mgrid.dy*M2.*expn2./(2i.*K);  % P03(z-mgrid.dy)
            F3(isnan(F3)) = 0;
            F3(real(Ktemp)<=0) = 0; 
    
            % 4     
            f  = ifft(ifftshift(F3)); % p0(z-mgrid.dy)  
            M3 = fftshift(fft(M_linear(:,I).*f));  % M(f1(z-mgrid.dy))
            F4 = F2.*expn2 + 0.25*mgrid.dy.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z-mgrid.dy)
            F4(isnan(F4)) = 0;
            F4(real(Ktemp)<=0) = 0;    
    
            % 5   
            f  = ifft(ifftshift(F4)); % p0(z-mgrid.dy)  
            M4 = fftshift(fft(M_linear(:,I).*f)); % M(f1(z-mgrid.dy))    
            F5 = excit_F.*exp(1i*K*mgrid.dy) +...
                 mgrid.dy/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*mgrid.dy)).*...
                 exp(1i*K*mgrid.dy)./(2i.*K); % P25(z-mgrid.dy)
            F5(isnan(F5)) = 0;
            F5(real(Ktemp)<=0) = 0; 
                
%% phase correction    
            wx2 = omegac.*omegac;
            e1 = exp(1i*(K - K2(:,I))*mgrid.dy);
            k_dif2 = (wx2./medium.c0./medium.c0 - wx2./c2./c2)*mgrid.dy;      
            k_dif  = (wx2./medium.c0./medium.c0 - wx2./c1./c1)*mgrid.dy; 
            clear  wx2  
            correction = excit_F.*exp(1i*K2(:,I)*mgrid.dy).*...
                        (1-e1.*(1+k_dif./4i./K)./(1-k_dif2./4i./K));                       
            F5 = F5 + correction; 
%             F5 = F5 ;     
            F5(isnan(F5)) = 0;
            F5(real(Ktemp)<=0) = 0;
            f = (ifft(ifftshift(F5)));
     
            T_rhoc = 2.*rho2.*c2./(rho1.*c1 + rho2.*c2); 
            f = f.*T_rhoc;
 
            ref_fundamental(:,I) = ref_fundamental(:,I) + f;  
                
%% reflection                
            p_ref2(:,I+1) = (ifft(ifftshift(excit_F))).*(T_rhoc-1);  
        end
    else
        p_ref = zeros(size(p_ref));
        f = 0;
        for I = 2:mgrid.num_y+1
            
            f = f + p_ref2(:,I-1); 
            excit_F = fftshift(fft(f));
    
            c1 = c_c(:,I-1);
            c2 = c_c(:,I);
            rho1 = rho_rho(:,I-1);
            rho2 = rho_rho(:,I);

            Ktemp1 = (omegac).^2./cw(:,I).^2 - evanescent_fiter*mgrid.kx.'.^2;
        
            % Simpson  
            % 1    
            M  = fftshift(fft(M_linear(:, I-1).*f));
            % P01(z+mgrid.dy/2)
            F1 = excit_F.*expn2 + 0.5*mgrid.dy*M.*expn2./(2i.*K); 
            F1(isnan(F1)) = 0;
            F1(real(Ktemp1)<=0) = 0; 
    
            % 2
            f  = ifft(ifftshift(F1));  % p0(z+mgrid.dy/2)
            M1 = fftshift(fft(M_linear(:, I-1).*f));
            % P12(z+mgrid.dy/2)
            F2 = excit_F.*expn2 + 0.25*mgrid.dy*expn2./(2i.*K).*(M + M1./expn2);
            clear M1  F1
            F2(isnan(F2)) = 0;
            F2(real(Ktemp1)<=0) = 0; 
    
            % 3   
            f  = ifft(ifftshift(F2)); % p1(z+mgrid.dy/2)  
            M2 = fftshift(fft(M_linear(:, I-1).*f));
            % P03(z+mgrid.dy)
            F3 = F2.*expn2 + 0.5*mgrid.dy*M2.*expn2./(2i.*K);  
            F3(isnan(F3)) = 0;
            F3(real(Ktemp1)<=0) = 0;   
    
            % 4     
            f  = ifft(ifftshift(F3)); % p0(z+mgrid.dy)  
            M3 = fftshift(fft(M_linear(:, I).*f));
            % P14(z+mgrid.dy)
            F4 = F2.*expn2 + 0.25*mgrid.dy.*(M2 + M3./expn2).*expn2./(2i.*K); 
            clear M3  F2
            F4(isnan(F4)) = 0;
            F4(real(Ktemp1)<=0) = 0;   
    
            % 5   
            f  = ifft(ifftshift(F4)); % p0(z+mgrid.dy)  
            M4 = fftshift(fft(M_linear(:, I).*f));    
            % P25(z+mgrid.dy)
            F5 = excit_F.*exp(1i*K*mgrid.dy) + ...
                 mgrid.dy/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*mgrid.dy)).*...
                 exp(1i*K*mgrid.dy)./(2i.*K); 
            clear M  M2  M4  F4
            F5(isnan(F5)) = 0;
            F5(real(Ktemp1)<=0) = 0; 
            %% phase correction 
            wx2 = omegac.*omegac;
            e1 = exp(1i*(K - K2(:,I))*mgrid.dy);
            k_dif2 = (wx2./medium.c0./medium.c0 - wx2./c2./c2)*mgrid.dy;      
            k_dif  = (wx2./medium.c0./medium.c0 - wx2./c1./c1)*mgrid.dy; 
            clear  wx2  
            correction = excit_F.*exp(1i*K2(:,I)*mgrid.dy).*...
                         (1-e1.*(1+k_dif./4i./K)./(1-k_dif2./4i./K));
         
            F5 = F5 + correction;  
    
            F5(isnan(F5)) = 0;
            F5(real(Ktemp1)<=0) = 0; 
            f = (ifft(ifftshift(F5)));
    
            T_rhoc = 2.*rho2.*c2./(rho1.*c1 + rho2.*c2);  %%layered medium 
            f = f.*T_rhoc;
    
            ref_fundamental(:,I) = ref_fundamental(:,I) + f;

            p_ref(:,I-1) = (ifft(ifftshift(excit_F))).*(T_rhoc-1);  
        end
        
        waitbar(iref/reflection_order)
    end
end
close (h)

end
