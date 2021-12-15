function [ref_fundamental2] = BReflection2D_fund(mgrid, medium, p_ref, ...
                             M_linear, K, cw, omega_c,...
                             reflection_order, excit_p, MDMS)

% DESCRIPTION:
% Simpson 2D integration in the forward direction at the fundamental
% frequency

% USAGE:
% [ref_fundamental2] = BReflection2D_fund(mgrid, medium, p_ref, ...
%                      M_linear, K, cw, omega_c,...
%                      reflection_order, excit_p, MDMS)

% INPUTS:
% mgrid             structure to define the computational domain
% medium            medium properties
% p_ref             stored reflected part
% M_linear          linear M_term
% K                 wave number along the y direction with a constant c
% cw                frequency-dependent speed of sound
% omega_c           transducer center frequency  
% reflection_order  the maximum order of reflection included in the simulation
% excit_ps          excitation signal
% MDMS              flag for reflection correction in backward propagation

% OUTPUTS:
% ref_fundamental2  reflections at the fundamental frequency

%%
% exponential term for the simpsion iteration
expn2 = exp(0.5*1i*K*mgrid.dy); 

% 1.05 is an empirical number, it can be higher if necessay
evanescent_fiter = 1.05;

% preallocate space for reflections
ref_fundamental = zeros(mgrid.num_x,mgrid.num_y+1); 
ref_fundamental2 = zeros(mgrid.num_x,mgrid.num_y+1); 

if length(cw)==1 
    cw = cw*ones(size(M_linear));
end

% c_c = medium.c;
% if length(medium.c) == 1
%     c_c = medium.c*ones(size(M_linear));
% end
    
rho_rho = medium.rho;
if length(medium.rho) == 1
    rho_rho = medium.rho*ones(size(M_linear));
end  
    
h = waitbar(0,'Reflection projection, please wait...');

for iref = 1:reflection_order
    
    if mod(iref,2)==1
        p_ref2 = zeros(mgrid.num_x,mgrid.num_y+1); 
        f = 0;
        for I = mgrid.num_y:-1:1 % the layer which the main waveform arrives 
                
                c1 = cw(:,I+1);
                c2 = cw(:,I);
                rho1 = rho_rho(:,I+1);
                rho2 = rho_rho(:,I);
    
                Ktemp = omega_c.^2./cw(:,I).^2 - evanescent_fiter*(mgrid.kx.').^2;    
                
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
                F5 = excit_F.*exp(1i*K*mgrid.dy) + ...
                     mgrid.dy/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*mgrid.dy)).*...
                     exp(1i*K*mgrid.dy)./(2i.*K); % P25(z-mgrid.dy)
                F5(isnan(F5)) = 0;
                F5(real(Ktemp)<=0) = 0;              
                f = (ifft(ifftshift(F5)));  
                ref_fundamental(:,I) = ref_fundamental(:,I) + f.*sqrt(rho2);  
                
%% reflection    
                T_rhoc = 2.*rho2.*c2./(rho1.*c1 + rho2.*c2); 
                p_ref2(:,I+1) = (ifft(ifftshift(excit_F))).*(T_rhoc-1);  
        end
        
    else
        p_ref = zeros(size(p_ref));
        f = 0;
        
        for I = 2:mgrid.num_y+1
                
                f = f + p_ref2(:,I-1); 
                excit_F = fftshift(fft(f));
    
                c1 = cw(:,I-1);
                c2 = cw(:,I);
                rho1 = rho_rho(:,I-1);
                rho2 = rho_rho(:,I);

                Ktemp1 = (omega_c).^2./cw(:,I).^2 - evanescent_fiter*mgrid.kx.'.^2;
        
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
                f = ifft(ifftshift(F5));
%                 ref_fundamental(:,I) = ref_fundamental(:,I) + f.*sqrt(rho2);
%% reflection
                T_rhoc = 2.*rho2.*c2./(rho1.*c1 + rho2.*c2);  %%layered medium  
                p_ref(:,I-1) = (ifft(ifftshift(excit_F))).*(T_rhoc-1);  
                
        end

        waitbar(iref/reflection_order)
    end
end
close (h)

% normalize the pressure field
f = excit_p.'./sqrt(rho_rho(:,end));
backward_dy = -mgrid.dy;
expn2 = exp(0.5*1i*K*backward_dy); 

if  MDMS == 1
    
    for I = mgrid.num_y:-1:1
        
        rho2 = rho_rho(:,I);
        c1 = cw(:,I+1);
        c2 = cw(:,I);
        rho1 = rho_rho(:,I+1);            
        T_rhoc = 2.*rho1.*c1./(rho1.*c1 + rho2.*c2);  %%layered medium 
       
        if  (max(max(abs(T_rhoc))) ~=1 || min(min(abs(T_rhoc))) ~=1) 
            f = f - p_ref2(:,I);       
        end
        excit_F = fftshift(fft(f));             

        Ktemp1 = (omega_c).^2./cw(:,I).^2 - evanescent_fiter*mgrid.kx.'.^2;
        
        % Simpson  
        % 1    
        M  = fftshift(fft(M_linear(:, I+1).*f));
        % P01(z+backward_dy/2)
        F1 = excit_F.*expn2 + 0.5*backward_dy*M.*expn2./(2i.*K); 
        F1(isnan(F1)) = 0;
        F1(real(Ktemp1)<=0) = 0; 
    
        % 2
        f  = ifft(ifftshift(F1));  % p0(z+backward_dy/2)
        M1 = fftshift(fft(M_linear(:, I+1).*f));
        % P12(z+backward_dy/2)
        F2 = excit_F.*expn2 + 0.25*backward_dy*expn2./(2i.*K).*(M + M1./expn2);
        clear M1  F1
        F2(isnan(F2)) = 0;
        F2(real(Ktemp1)<=0) = 0;     
        % 3   
        f  = ifft(ifftshift(F2)); % p1(z+backward_dy/2)  
        M2 = fftshift(fft(M_linear(:, I+1).*f));
        % P03(z+backward_dy)
        F3 = F2.*expn2 + 0.5*backward_dy*M2.*expn2./(2i.*K);  
        F3(isnan(F3)) = 0;
        F3(real(Ktemp1)<=0) = 0;     
        % 4     
        f  = ifft(ifftshift(F3)); % p0(z+backward_dy)  
        M3 = fftshift(fft(M_linear(:, I).*f));
        % P14(z+backward_dy)
        F4 = F2.*expn2 + 0.25*backward_dy.*(M2 + M3./expn2).*expn2./(2i.*K); 
        clear M3  F2
        F4(isnan(F4)) = 0;
        F4(real(Ktemp1)<=0) = 0;     
        % 5   
        f  = ifft(ifftshift(F4)); % p0(z+backward_dy)  
        M4 = fftshift(fft(M_linear(:, I).*f));    
        % P25(z+backward_dy)
        F5 = excit_F.*exp(1i*K*backward_dy) +...
             backward_dy/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*backward_dy)).*exp(1i*K*backward_dy)./(2i.*K); 
        clear M  M2  M4  F4
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
        f = ifft(ifftshift(F5));
%         excit_F = F5;
        clear F5
        
        % recover pressure from the normalized wave field
        ref_fundamental2(:,I) = ...
            ref_fundamental(:,I) + f.*sqrt(rho2);
    end
    ref_fundamental2(:,end) = excit_p; 
    
else       
        
    ref_fundamental2 = ref_fundamental;    
    
end


end
