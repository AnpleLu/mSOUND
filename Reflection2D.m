function [p_reflection] = Reflection2D(mgrid, medium, p_ref, ...
                          M_linear, K, cw, sensor_mask,...
                          reflection_order)
% DESCRIPTION:
% Simpson-like intergration for the reflection projection

% USAGE:
% [p_reflection] = Reflection2D(mgrid, medium, p_ref, ...
%                  M_linear, K, cw, sensor_mask, reflection_order)

% INPUTS:
% mgrid             structure to define the computational domain
% medium            medium properties
% p_ref             stored reflected part
% M_linear          linear part of the M term on RHS of solution equation 
%                   obtained with function Mterm1D 
% K                 wave number along the y direction with constant c
% cw                frequency-dependent speed of sound
% sensor_mask       a set of Cartesian points where the pressure is recorded
% reflection_order  the maximum order of reflection included in the simulation

% OUTPUTS:
% p_reflection      reflection     
%%
% precalculated const. for the iteration         
expn2 = exp(0.5*1i*K*mgrid.dy); 
                
% allocate space for the reflections
p_reflection = zeros(mgrid.num_t, sum(sum(sensor_mask)));

rho_rho = medium.rho;
if length(medium.rho) == 1
    rho_rho = medium.rho*ones(mgrid.num_x, mgrid.num_y+1);
end 

% waitbar
h = waitbar(0,'Reflection projection, please wait...');

% 1.05 is an empirical number, it can be higher if necessay
evanescent_fiter = 1.05;

% Simpson iteration
for iref = 1:reflection_order
    
    if mod(iref,2)==1
        
        % initial an array to store reflected parts
        p_ref2 = zeros(mgrid.num_t,mgrid.num_x,mgrid.num_y+1); 
        f = 0;
        for I = mgrid.num_y:-1:1 % the layer at which the main waveform arrives
          
            % reflection for propagation
            f = f + p_ref(:,:,I+1);   
            
            % low pass filter
            Ktemp = mgrid.w'.^2*ones(1,mgrid.num_x)./cw(:,:,I).^2 -...
                    evanescent_fiter*ones(mgrid.num_t,1)*mgrid.kx.^2;    
            
                
            excit_F = fftshift(fft2(f));
            % Simpson   
            % 1    
            Ft  = fftshift(fft(f,[],1),1);
            M = fftshift(fft((M_linear(:,:,I+1)).*Ft,[],2),2);           % M(f(z))
            clear Ft  f
            F1 = excit_F.*expn2 + 0.5*mgrid.dy*M.*expn2./(2i.*K);  % P01(z-mgrid.dy/2)
            F1(isnan(F1)) = 0;
            F1(real(Ktemp)<=0) = 0; 

            % 2
            f   = real(ifft2(ifftshift(F1)));  % p0(z-mgrid.dy/2)
            Ft  = fftshift(fft(f,[],1),1);     % fft(f) 
            M1 = fftshift(fft((M_linear(:,:,I+1)).*Ft,[],2),2);           % M(f(z))
            F2 = excit_F.*expn2 + 0.25*mgrid.dy*expn2./(2i.*K).*(M + M1./expn2);% P12(z-mgrid.dy/2)
            F2(isnan(F2)) = 0;
            F2(real(Ktemp)<=0) = 0; 
    
            % 3   
            f   = real(ifft2(ifftshift(F2))); % p1(z-mgrid.dy/2)  
            Ft  = fftshift(fft(f,[],1),1);    % fft(f1) 
            M2  = fftshift(fft((M_linear(:,:,I+1)).*Ft,[],2),2);     % M(f1(z-mgrid.dy/2))
            F3 = F2.*expn2 + 0.5*mgrid.dy*M2.*expn2./(2i.*K);  % P03(z-mgrid.dy)
            F3(isnan(F3)) = 0;
            F3(real(Ktemp)<=0) = 0; 
    
            % 4     
            f   = real(ifft2(ifftshift(F3))); % p0(z-mgrid.dy)  
            Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
            M3  = fftshift(fft((M_linear(:,:,I)).*Ft,[],2),2);  % M(f1(z-mgrid.dy))
            F4 = F2.*expn2 + 0.25*mgrid.dy.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z-mgrid.dy)
            F4(isnan(F4)) = 0;
            F4(real(Ktemp)<=0) = 0;    
      
            % 5   
            f   = real(ifft2(ifftshift(F4))); % p0(z-mgrid.dy)  
            Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
            M4  = fftshift(fft((M_linear(:,:,I)).*Ft,[],2),2); % M(f1(z-mgrid.dy))    
            F5 = excit_F.*exp(1i*K*mgrid.dy) +...
                 mgrid.dy/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*mgrid.dy)).*...
                 exp(1i*K*mgrid.dy)./(2i.*K); % P25(z-mgrid.dy)
            F5(isnan(F5)) = 0;
            F5(real(Ktemp)<=0) = 0; 
            f = real(ifft2(ifftshift(F5))); 
           
            % recover from the normalized wave field
            rho2 = (repmat(rho_rho(:,I)', mgrid.num_t,1)); 
            p_I = sqrt(rho2).*f;    
    
            if sum(sum(sensor_mask(:,I)))~=0
                % calculate the starting position of recorded signal
                sensor_mask_I1 = sum(sum(sensor_mask(:,1:I-1)))+1;
       
                % calculate the ending position of recorded signal
                sensor_mask_I2 = sensor_mask_I1 -1 + ...
                sum(sum(sensor_mask(:,I)));
    
                % save the time-domain signals
                p_reflection(:,sensor_mask_I1:sensor_mask_I2) = ...
                   p_reflection(:,sensor_mask_I1:sensor_mask_I2) +...
                   p_I(:,sensor_mask(:,I)~=0);
            
            end            

            rho1 = (repmat(rho_rho(:,I+1)', mgrid.num_t,1));  
            % transmission coefficient
            T_rhoc = 2.*rho2.*cw(:,:,I)./(rho1.*cw(:,:,I+1) + rho2.*cw(:,:,I));  
            clear  rho1 rho2
            % store the reflected part
            p_ref2(:,:,I+1) = real(ifft2(ifftshift(excit_F))).*(T_rhoc-1);          
            clear T_rhoc      
        end
        
    else
        % initialize an array to store the reflected part
        p_ref = zeros(size(p_ref));
        f = 0;
        for I = 2:mgrid.num_y+1 % the layer at which the main waveform arrives

            Ktemp = mgrid.w'.^2*ones(1,mgrid.num_x)./cw(:,:,I).^2 -...
                    evanescent_fiter*ones(mgrid.num_t,1)*mgrid.kx.^2;    
     
            f = f + squeeze(p_ref2(:,:,I-1));     
            M_tot = M_linear(:,:,I-1);
            excit_F = fftshift(fft2(f));
 
            % 1    
            Ft = fftshift(fft(f,[],1),1);
            M = fftshift(fft((M_tot).*Ft,[],2),2);  % M(f(z))
            F1 = excit_F.*expn2 + 0.5*mgrid.dy*M.*expn2./(2i.*K);  % P01(z-mgrid.dy/2)
            F1(isnan(F1)) = 0;
            F1(real(Ktemp)<=0) = 0; 

            % 2
            f  = real(ifft2(ifftshift(F1)));  % p0(z-mgrid.dy/2)
            Ft = fftshift(fft(f,[],1),1);     % fft(f) 
            M1 = fftshift(fft((M_tot).*Ft,[],2),2);   % M(f(z))
            F2 = excit_F.*expn2 + 0.25*mgrid.dy*expn2./(2i.*K).*(M + M1./expn2);% P12(z-mgrid.dy/2)
            F2(isnan(F2)) = 0;
            F2(real(Ktemp)<=0) = 0; 
    
            % 3   
            f  = real(ifft2(ifftshift(F2))); % p1(z-mgrid.dy/2)  
            Ft = fftshift(fft(f,[],1),1);    % fft(f1) 
            M2 = fftshift(fft((M_tot).*Ft,[],2),2);     % M(f1(z-mgrid.dy/2))
            F3 = F2.*expn2 + 0.5*mgrid.dy*M2.*expn2./(2i.*K);  % P03(z-mgrid.dy)
            F3(isnan(F3)) = 0;
            F3(real(Ktemp)<=0) = 0; 
    
            % 4  
            M_tot = M_linear(:,:,I);
            f  = real(ifft2(ifftshift(F3))); % p0(z-mgrid.dy)  
            Ft = fftshift(fft(f,[],1),1);    % fft(f2) 
            M3 = fftshift(fft((M_tot).*Ft,[],2),2);  % M(f1(z-mgrid.dy))
            F4 = F2.*expn2 + 0.25*mgrid.dy.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z-mgrid.dy)
            F4(isnan(F4)) = 0;
            F4(real(Ktemp)<=0) = 0;    
    
            % 5   
            f  = real(ifft2(ifftshift(F4))); % p0(z-mgrid.dy)  
            Ft = fftshift(fft(f,[],1),1);    % fft(f2) 
            M4 = fftshift(fft((M_tot).*Ft,[],2),2); % M(f1(z-mgrid.dy))    
            F5 = excit_F.*exp(1i*K*mgrid.dy) +...
                 mgrid.dy/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*mgrid.dy)).*...
                 exp(1i*K*mgrid.dy)./(2i.*K); % P25(z-mgrid.dy)
            F5(isnan(F5)) = 0;
            F5(real(Ktemp)<=0) = 0; 
     
            f = real(ifft2(ifftshift(F5)));     
            rho2 = (repmat(rho_rho(:,I)', mgrid.num_t,1));
            p_I = sqrt(rho2).*f;   
            
            if sum(sum(sensor_mask(:,I)))~=0
                % calculate the starting position of recorded signal
                sensor_mask_I1 = sum(sum(sensor_mask(:,1:I-1)))+1;
       
                % calculate the ending position of recorded signal
                sensor_mask_I2 = sensor_mask_I1 -1 + ...
                sum(sum(sensor_mask(:,I)));
    
                 % save the time-domain signals
                 p_reflection(:,sensor_mask_I1:sensor_mask_I2) = ...
                     p_reflection(:,sensor_mask_I1:sensor_mask_I2) +...
                     p_I(:,sensor_mask(:,I)~=0);
            
            end
            
            c1 = cw(:,:,I-1);
            c2 = cw(:,:,I);    
            rho1 = (repmat(rho_rho(:,I-1)', mgrid.num_t,1));
            T_rhoc = 2.*rho2.*c2./(rho1.*c1 + rho2.*c2);  %%layered medium
            clear c1 c2 rho1 
            
            % store the reflected parts
            p_ref(:,:,I-1) = real(ifft2(ifftshift(excit_F))).*(T_rhoc-1);
            clear T_rhoc
        end        
         
    end
    waitbar(iref/reflection_order)
end
% close waitbar
close(h)

end


