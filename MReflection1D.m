function [p_reflection] = MReflection1D(mgrid, medium, p_ref, ...
                          M_linear, K, K2, cw, sensor_mask,...
                          reflection_order)
  

% DESCRIPTION:
% Simpson-like intergration for reflection propagation with tranmission corrections

% USAGE:
% [p_reflection] = MReflection1D(mgrid, medium, p_ref, ...
%                  M_linear, K, K2, cw, sensor_mask, ...
%                  reflection_order)

% INPUTS:
% mgrid              structure to define the computational domain
% medium             medium properties
% p_ref              stored reflected part
% M_linear           linear part of the M term on RHS of solution equation 
%                    obtained with function Mterm1D_MMDM            
% K                  wavenumber in the z direction with constant c
% K2                 wavenumber in the z direction with heterogeneous c
% cw                 frequency-dependent speed of sound
% sensor_mask        a set of Cartesian points where the pressure is recorded
% reflection_order   the maximum order of reflection included in the simulation

% OUTPUTS:
% p_reflection       reflection          
             
%%  

% exponential term for the forward propagation
expn2 = exp(0.5*1i*K*mgrid.dx); 

rho_rho = medium.rho;
if length(medium.rho) == 1
    rho_rho = medium.rho*ones(mgrid.num_x+1);
end  

% preallocate space for the reflections
p_reflection = zeros(mgrid.num_t, sum(sum(sensor_mask)));

h = waitbar(0,'reflection projection, please wait...');

for iref = 1:reflection_order
    
    if mod(iref,2)==1
        p_ref2 = zeros(mgrid.num_t,mgrid.num_x+1); 
        f = 0;
    for I = mgrid.num_x:-1:1     
        % Simpson integral algorithm
        c1 = cw(:,I+1);
        c2 = cw(:,I);    
        rho1 = rho_rho(I+1);
        rho2 = rho_rho(I);
 
        % 1 
        f = f + p_ref(:,I+1);  
        excit_F = fftshift(fft(f));
        Ft  = fftshift(fft(f,[],1),1);
        % M(f(z))
        M = M_linear(:,I+1).*Ft;
        F1 = excit_F.*expn2 + 0.5*mgrid.dx*M.*expn2./(2i.*K);  % P01(z+mgrid.dx/2)
        F1(isnan(F1)) = 0;

        % 2
        f = real(ifft(ifftshift(F1)));  % p0(z+mgrid.dx/2) 
        clear F1
        Ft  = fftshift(fft(f,[],1),1);     % fft(f) 
        % M(f0(z+mgrid.dx/2))  M1
        M1 = M_linear(:,I+1).*Ft;
        F2 = excit_F.*expn2 + 0.25*mgrid.dx*expn2./(2i.*K).*(M + M1./expn2);% P12(z+mgrid.dx/2)
        F2(isnan(F2)) = 0;
  
        % 3   
        f   = real(ifft(ifftshift(F2))); % p1(z+mgrid.dx/2)   
        Ft  = fftshift(fft(f,[],1),1);    % fft(f1) 
        % M(f1(z+mgrid.dx/2)) M2
        M2 = M_linear(:,I+1).*Ft;
        F3 = F2.*expn2 + 0.5*mgrid.dx*M2.*expn2./(2i.*K);  % P03(z+mgrid.dx)
        F3(isnan(F3)) = 0;
    
        % 4     
        f   = real(ifft(ifftshift(F3))); % p0(z+mgrid.dx)  
        clear F3
        Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
        % M(f0(z+mgrid.dx))  M3
        M3 = M_linear(:,I).*Ft;
        F4 = F2.*expn2 + 0.25*mgrid.dx.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z+mgrid.dx)
        clear F2
        F4(isnan(F4)) = 0;  
      
        % 5   
        f   = real(ifft(ifftshift(F4))); % p0(z+mgrid.dx)
        clear F4
        Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
        % M(f1(z+mgrid.dx)) M4
        M4 = M_linear(:,I).*Ft;
        F5 = excit_F.*exp(1i*K*mgrid.dx) +...
             mgrid.dx/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*mgrid.dx)).*...
             exp(1i*K*mgrid.dx)./(2i.*K); % P25(z+mgrid.dx)
     
        clear M  M2  M4  
        F5(isnan(F5)) = 0;

        %% phase correction
        wx2 = (mgrid.w.*mgrid.w).'; 
        e1 = exp(1i*(K - K2(:,I))*mgrid.dx);

        k_dif2 = (wx2./medium.c0./medium.c0 - wx2./c2./c2)*mgrid.dx;      
        k_dif  = (wx2./medium.c0./medium.c0 - wx2./c1./c1)*mgrid.dx; 
    
        clear  wx2  
        correction = excit_F.*exp(1i*K2(:,I)*mgrid.dx).*...
                     (1-e1.*(1+k_dif./4i./K)./(1-k_dif2./4i./K));
                
        F5 = F5 + correction;  
        clear correction k_dif2  k_dif  el

        F5(isnan(F5)) = 0;
   
        %% amplitude correction
        T_rhoc = 2.*rho2.*c2./(rho1.*c1 + rho2.*c2);  %%layered medium
        f = real(ifft(ifftshift(F5)));
        f = f.*T_rhoc;
        F5 = fftshift(fft2(f));  

        f = real(ifft(ifftshift(F5))); % p0(z+dz) 
        clear F5
        %% save data
    
        if sum(sensor_mask(I))~=0
            
            % calculate the starting position of recorded signal
            sensor_mask_I1 = sum((sensor_mask(1:I-1)))+1;
        
            % calculate the ending position of recorded signal
            sensor_mask_I2 = sensor_mask_I1 -1 + ...
            sum(sensor_mask(I));
        
            % save the time-domain signals
            p_reflection(:,sensor_mask_I1:sensor_mask_I2) = ...
                p_reflection(:,sensor_mask_I1:sensor_mask_I2) + f;
        
        end   
%% reflection        
        p_ref2(:,I+1) = real(ifft(ifftshift(excit_F))).*(T_rhoc-1);         
        
    
    end
        
    else
        p_ref = zeros(size(p_ref));
        f = 0;
    for I = 2:mgrid.num_x+1  
        
         % Simpson integral algorithm
        c1 = cw(:,I-1);
        c2 = cw(:,I);    
        rho1 = rho_rho(I-1);
        rho2 = rho_rho(I);
        
        %Simpson integral algorithm

        % 1   
        f = f + squeeze(p_ref2(:,I-1)); 
        excit_F = fftshift(fft(f));
        Ft  = fftshift(fft(f,[],1),1);
        % M(f(z))
        M = M_linear(:,I-1).*Ft ;
        F1 = excit_F.*expn2 + 0.5*mgrid.dx*M.*expn2./(2i.*K);  % P01(z+mgrid.dx/2)
        F1(isnan(F1)) = 0;
    
        % 2
        f = real(ifft(ifftshift(F1)));  % p0(z+mgrid.dx/2) 
        clear F1
        Ft  = fftshift(fft(f,[],1),1);     % fft(f) 
        % M(f0(z+mgrid.dx/2))  M1
        M1 = M_linear(:,I-1).*Ft;
        F2 = excit_F.*expn2 + 0.25*mgrid.dx*expn2./(2i.*K).*(M + M1./expn2);% P12(z+mgrid.dx/2)
        F2(isnan(F2)) = 0;
   
        % 3   
        f   = real(ifft(ifftshift(F2))); % p1(z+mgrid.dx/2)   
        Ft  = fftshift(fft(f,[],1),1);    % fft(f1) 
        % M(f1(z+mgrid.dx/2)) M2
        M2 = M_linear(:,I-1).*Ft;
        F3 = F2.*expn2 + 0.5*mgrid.dx*M2.*expn2./(2i.*K);  % P03(z+mgrid.dx)
        F3(isnan(F3)) = 0;
    
        % 4     
        f   = real(ifft(ifftshift(F3))); % p0(z+mgrid.dx)  
        clear F3
        Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
        % M(f0(z+mgrid.dx))  M3
        M3 = M_linear(:,I).*Ft;
        F4 = F2.*expn2 + 0.25*mgrid.dx.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z+mgrid.dx)
        clear F2
        F4(isnan(F4)) = 0;  
      
        % 5   
        f   = real(ifft(ifftshift(F4))); % p0(z+mgrid.dx)
        clear F4
        Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
        % M(f1(z+mgrid.dx)) M4
        M4 = M_linear(:,I).*Ft;
        F5 = excit_F.*exp(1i*K*mgrid.dx) +...
             mgrid.dx/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*mgrid.dx)).*...
             exp(1i*K*mgrid.dx)./(2i.*K); % P25(z+mgrid.dx)
     
        clear M  M2  M4  
        F5(isnan(F5)) = 0;
   
        %% phase correction
        wx2 = (mgrid.w.*mgrid.w).'; 
        e1 = exp(1i*(K - K2(:,I))*mgrid.dx);

        k_dif2 = (wx2./medium.c0./medium.c0 - wx2./c2./c2)*mgrid.dx;      
        k_dif  = (wx2./medium.c0./medium.c0 - wx2./c1./c1)*mgrid.dx; 
    
        clear  wx2  
        correction = excit_F.*exp(1i*K2(:,I)*mgrid.dx).*...
                     (1-e1.*(1+k_dif./4i./K)./(1-k_dif2./4i./K));
                
        F5 = F5 + correction;  
        clear correction k_dif2  k_dif  el

        F5(isnan(F5)) = 0;
   
        %% amplitude correction
        T_rhoc = 2.*rho2.*c2./(rho1.*c1 + rho2.*c2);  %%layered medium
        f = real(ifft(ifftshift(F5)));
        f = f.*T_rhoc;
        F5 = fftshift(fft2(f));          
  
 %%       
        f = real(ifft(ifftshift(F5))); % p0(z+mgrid.dx) 
        clear F5
%% save data
        if sum(sensor_mask(I))~=0
            
            % calculate the starting position of recorded signal
            sensor_mask_I1 = sum((sensor_mask(1:I-1)))+1;
        
            % calculate the ending position of recorded signal
            sensor_mask_I2 = sensor_mask_I1 -1 + ...
            sum(sensor_mask(I));
        
            % save the time-domain signals
            p_reflection(:,sensor_mask_I1:sensor_mask_I2) = ...  
                p_reflection(:,sensor_mask_I1:sensor_mask_I2) + f;
        end

     
        clear p_I
        p_ref(:,I-1) = real(ifft(ifftshift(excit_F))).*(T_rhoc-1); 
    end 
         
    end
    waitbar(iref/reflection_order)
end
close(h)

end


