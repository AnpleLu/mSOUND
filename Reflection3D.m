function [p_reflection] = Reflection3D(mgrid, medium, p_ref, ...
                          M_linear, K, cw, sensor_mask,...
                          reflection_order)
% DESCRIPTION:
% Simpson-like integration for the reflection projection with correction

% USAGE:
% [p_reflection] = Reflection3D(mgrid, medium, p_ref, ...
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

% exponential term for the forward propagation
expn2 = exp(0.5*1i*K*mgrid.dz);  

rho_rho = medium.rho;
if length(medium.rho) == 1
    rho_rho = medium.rho*ones(mgrid.num_x, mgrid.num_y, mgrid.num_z+1);
end   

% preallocation
rho1 = zeros(1, mgrid.num_x, mgrid.num_y);
rho2 = zeros(1, mgrid.num_x, mgrid.num_y);   

% 1.05 is an empirical number, it can be higher if necessay
evanescent_fiter = 1.05;

%  prepare for the calculation of evanescent wave
Ktemp1_kx = zeros(1, mgrid.num_x);
Ktemp1_ky = zeros(1,1,mgrid.num_y);
Ktemp1_kx(1, :)  = squeeze(mgrid.kx);
Ktemp1_ky(1,1,:) = squeeze(mgrid.ky);
    
% allocate space for the reflections
p_reflection = zeros(mgrid.num_t, sum(sum(sum(sensor_mask))));

% waitar
h = waitbar(0,'Reflection projection, please wait...');

for iref = 1:reflection_order  
    
    if mod(iref,2)==1
        p_ref2 = zeros(mgrid.num_t,mgrid.num_x,mgrid.num_y, mgrid.num_z+1); 
        f = 0;
        for I = mgrid.num_z:-1:1 % the layer at which the main waveform arrives
            
            % filter for filtering out evanescent wave
            Ktemp1 = (repmat((mgrid.w.'),1,mgrid.num_x, mgrid.num_y)./...
                     cw(:,:,:,I)).^2 - ...
                     evanescent_fiter*repmat(Ktemp1_kx.^2, mgrid.num_t, 1, mgrid.num_y) -...
                     evanescent_fiter*repmat(Ktemp1_ky.^2, mgrid.num_t, mgrid.num_x, 1);
       
            % 1 
            f = f + p_ref(:,:,:,I+1);   
            excit_F = fftshift(fftn(f));
            Ft  = fftshift(fft(f,[],1),1);
            % M(f(z))
            M = fftshift(fft(fftshift(fft((M_linear(:,:,:,I+1)).*Ft,[],2),2),[],3),3);
            clear   Ft  F2t
            F1 = excit_F.*expn2 + 0.5*mgrid.dz*M.*expn2./(2i.*K);  % P01(z+mgrid.dz/2)
            F1(isnan(F1)) = 0;
            F1(real(Ktemp1)<=0) = 0; 
   
            % 2
            f = real(ifftn(ifftshift(F1)));  % p0(z+mgrid.dz/2) 
            clear F1
            Ft  = fftshift(fft(f,[],1),1);     % fft(f) 
            % M(f0(z+mgrid.dz/2))  M1
            M1 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I+1)).*Ft,[],2),2),[],3),3);
            clear   Ft  F2t
            F2 = excit_F.*expn2 + 0.25*mgrid.dz*expn2./(2i.*K).*(M + M1./expn2);% P12(z+mgrid.dz/2)
            clear M1
            F2(isnan(F2)) = 0;
            F2(real(Ktemp1)<=0) = 0;   
  
            % 3   
            f   = real(ifftn(ifftshift(F2))); % p1(z+mgrid.dz/2)   
            Ft  = fftshift(fft(f,[],1),1);    % fft(f1) 
            % M(f1(z+mgrid.dz/2)) M2
            M2 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I+1)).*Ft,[],2),2),[],3),3);
            clear   Ft  F2t
            F3 = F2.*expn2 + 0.5*mgrid.dz*M2.*expn2./(2i.*K);  % P03(z+mgrid.dz)
            F3(isnan(F3)) = 0;
            F3(real(Ktemp1)<=0) = 0;  
    
            % 4     
            f   = real(ifftn(ifftshift(F3))); % p0(z+mgrid.dz)   
            Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
            % M(f0(z+mgrid.dz))  M3
            M3 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I)).*Ft,[],2),2),[],3),3);
            F4 = F2.*expn2 + 0.25*mgrid.dz.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z+mgrid.dz)
            clear M3   Ft  F2t
            F4(isnan(F4)) = 0;
            F4(real(Ktemp1)<=0) = 0;    
      
            % 5   
            f   = real(ifftn(ifftshift(F4))); % p0(z+mgrid.dz) 
            clear F4
            Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
            % M(f1(z+mgrid.dz)) M4
            M4 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I)).*Ft,[],2),2),[],3),3);
            F5 = excit_F.*exp(1i*K*mgrid.dz) + ...
                 mgrid.dz/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*mgrid.dz)).*...
                 exp(1i*K*mgrid.dz)./(2i.*K); % P25(z+mgrid.dz)
     
            clear M   M2   M4  Ft  F2t
            F5(isnan(F5)) = 0;
            F5(real(Ktemp1)<=0) = 0; 
            f = real(ifftn(ifftshift(F5)));
            clear F5
            
            % recover the normalized wave field
            rho2(1,:,:) = squeeze(rho_rho(:,:, I)); 
            p_I = sqrt(repmat(rho2, mgrid.num_t, 1,1)).*f;    
    
            if sum(sum(sum(sensor_mask(:,:,I))))~=0
                % calculate the starting position of recorded signal
                sensor_mask_I1 = sum(sum(sum(sensor_mask(:,:,1:I-1)))) + 1;
               
                % calculate the ending position of recorded signal
                sensor_mask_I2 = sensor_mask_I1 -1 +  ...
                sum(sum(sum(sensor_mask(:,:,I))));
    
                % save the time-domain signals
                p_reflection(:,sensor_mask_I1:sensor_mask_I2) = ...
                    p_reflection(:,sensor_mask_I1:sensor_mask_I2) +...
                    p_I(:,sensor_mask(:,:,I)~=0);
            end
            clear p_I
           
            rho1(1,:,:) = rho_rho(:,:,I+1);
            rho2(1,:,:) = rho_rho(:,:, I);        
            rho11 = repmat(rho1, mgrid.num_t, 1,1); 
            rho22 = repmat(rho2, mgrid.num_t, 1,1);            
            T_rhoc = 2.*rho22.*cw(:,:,:,I)./...
                     (rho11.*cw(:,:,:,I+1) + rho22.*cw(:,:,:,I));          
            clear rho11 rho22 rho1 
            
            p_ref2(:,:,:, I+1) = real(ifftn(ifftshift(excit_F))).*(T_rhoc-1);
                      
            clear T_rhoc           
        end

    else
        p_ref = zeros(size(p_ref));
        f = 0;
        for I = 2:mgrid.num_z+1 % the layer at which the main waveform arrives

            % used for filtering out evanescent wave
            Ktemp1 = (repmat((mgrid.w.'),1,mgrid.num_x, mgrid.num_y)./...
                      cw(:,:,:,I)).^2 - ...
                      evanescent_fiter*repmat(Ktemp1_kx.^2, mgrid.num_t, 1, mgrid.num_y)-...
                      evanescent_fiter*repmat(Ktemp1_ky.^2, mgrid.num_t, mgrid.num_x, 1);

            % 1        
            f = f + squeeze(p_ref2(:,:,:,I-1));         
            excit_F = fftshift(fftn(f));

            Ft  = fftshift(fft(f,[],1),1);
            % M(f(z))
            M = fftshift(fft(fftshift(fft((M_linear(:,:,:,I-1)).*Ft,[],2),2),[],3),3);
            clear   Ft  F2t
            F1 = excit_F.*expn2 + 0.5*mgrid.dz*M.*expn2./(2i.*K);  % P01(z+mgrid.dz/2)
            F1(isnan(F1)) = 0;
            F1(real(Ktemp1)<=0) = 0; 
   
            % 2
            f = real(ifftn(ifftshift(F1)));  % p0(z+mgrid.dz/2) 
            clear F1
            Ft  = fftshift(fft(f,[],1),1);     % fft(f) 
            % M(f0(z+mgrid.dz/2))  M1
            M1 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I-1)).*Ft,[],2),2),[],3),3);
            clear   Ft  F2t
            F2 = excit_F.*expn2 + 0.25*mgrid.dz*expn2./(2i.*K).*(M + M1./expn2);% P12(z+mgrid.dz/2)
            clear M1
            F2(isnan(F2)) = 0;
            F2(real(Ktemp1)<=0) = 0;   
  
            % 3   
            f   = real(ifftn(ifftshift(F2))); % p1(z+mgrid.dz/2)   
            Ft  = fftshift(fft(f,[],1),1);    % fft(f1) 
            % M(f1(z+mgrid.dz/2)) M2
            M2 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I-1)).*Ft,[],2),2),[],3),3);
            clear   Ft  F2t
            F3 = F2.*expn2 + 0.5*mgrid.dz*M2.*expn2./(2i.*K);  % P03(z+mgrid.dz)
            F3(isnan(F3)) = 0;
            F3(real(Ktemp1)<=0) = 0;  
    
            % 4     
            f   = real(ifftn(ifftshift(F3))); % p0(z+mgrid.dz)   
            Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
            % M(f0(z+mgrid.dz))  M3
            M3 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I)).*Ft,[],2),2),[],3),3);
            F4 = F2.*expn2 + 0.25*mgrid.dz.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z+mgrid.dz)
            clear M3   Ft  F2t
            F4(isnan(F4)) = 0;
            F4(real(Ktemp1)<=0) = 0;    
      
            % 5   
            f   = real(ifftn(ifftshift(F4))); % p0(z+mgrid.dz) 
            clear F4
            Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
            % M(f1(z+mgrid.dz)) M4
            M4 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I)).*Ft,[],2),2),[],3),3);
            F5 = excit_F.*exp(1i*K*mgrid.dz) + ...
                 mgrid.dz/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*mgrid.dz)).*...
                 exp(1i*K*mgrid.dz)./(2i.*K); % P25(z+mgrid.dz)
     
            clear M   M2   M4  Ft  F2t
            F5(isnan(F5)) = 0;
            F5(real(Ktemp1)<=0) = 0; 
            f = real(ifftn(ifftshift(F5)));
            clear F5
            
            % recover the normalized wave field
            rho2(1,:,:) = squeeze(rho_rho(:,:, I));  
            p_I = sqrt(repmat(rho2, mgrid.num_t, 1,1)).*f;    
    
            if sum(sum(sum(sensor_mask(:,:,I))))~=0
               % calculate the starting position of recorded signal
               sensor_mask_I1 = sum(sum(sum(sensor_mask(:,:,1:I-1)))) + 1;
               
               % calculate the ending position of recorded signal
               sensor_mask_I2 = sensor_mask_I1 -1 +  ...
               sum(sum(sum(sensor_mask(:,:,I))));
    
               % save the time-domain signals
               p_reflection(:,sensor_mask_I1:sensor_mask_I2) = ...
                    p_reflection(:,sensor_mask_I1:sensor_mask_I2) +...
                    p_I(:,sensor_mask(:,:,I)~=0);
            end
           clear p_I         
           
           rho1(1,:,:) = rho_rho(:,:,I-1);
           rho2(1,:,:) = rho_rho(:,:, I);        
           rho11 = repmat(rho1, mgrid.num_t, 1,1); 
           rho22 = repmat(rho2, mgrid.num_t, 1,1); 
           % transmission coefficient
           T_rhoc = 2.*rho22.*cw(:,:,:,I)./...
                   (rho11.*cw(:,:,:,I-1) + rho22.*cw(:,:,:,I));          
           clear rho11 rho22 rho1 
           p_ref(:,:,:, I-1) = real(ifftn(ifftshift(excit_F))).*(T_rhoc-1);
           clear T_rhoc         
        end

    end
    waitbar(iref/reflection_order)
end
close(h) % close waitbar

end


