
function [P_second] = Forward3D_sec(mgrid, medium, P_fundamental, omega_c,...
                                    varargin)
% DESCRIPTION:
% Simpson-like 3D intergration in the forward direction at the second-harmonic
% frequency

% USAGE:
% [P_second] = Forward3D_sec(mgrid, medium, P_fundamental, omega_c,...
%                            varargin)

% INPUTS:
% mgrid               input structure to define the computational domain
% medium              medium properties
% P_fundamental       fundamental pressure dirtribution through the domain
% omega_c              center frequency 2*pi*fc

% OPTIONAL INOUTS:
% 'NRL'              add non-reflecting layer at the boundary to reduce the
% spatial aliasing error
% 'correction'       add correction due to strong heterogeneities


% OUTPUTS:
% P_second          second-harmonic pressure dirtribution through the domain

%%
MDMC = 0;  % flag for correction
MDMA = 0;  % flag for non-reflecting layer
if ~isempty(varargin)
    for input_index = 1:length(varargin)
        switch varargin{input_index}
            case 'correction'
            MDMC = 1;    % flag for transmission correction
            case 'NRL'
            MDMA = 1;    % flag for absorption
        end
    end
end

%% non-reflecting layer
if MDMA ==1
    M_gamma = repmat(1./(cosh(medium.NRL_alpha.*mgrid.abx_vec.')).^2, 1, mgrid.num_y) + ...
              repmat(1./(cosh(medium.NRL_alpha.*mgrid.aby_vec)).^2, mgrid.num_x, 1);
    M_gamma = M_gamma.*medium.NRL_gamma.*1i*omega_c*2;
else 
    M_gamma = 0;
end 

% 1.05 is an empirical number, it can be higher if necessay
evanescent_fiter = 1.05;

if MDMC == 0
    % calculate the M term on RHS of solution equation
    [M_linear, Ktemp, cw2, cw] = Mterm3D_sec(mgrid, medium, omega_c);
    
    % add non-reflecting layer
    M_linear = M_linear + repmat(M_gamma,1,1,mgrid.num_z+1);
    
    % read paper evaluation of a wave-vector...when omega is positive, K is negative.
    kxy  = squeeze((mgrid.kx.'*ones(1, mgrid.num_y)).^2 +...
           (ones(mgrid.num_x, 1)*mgrid.ky).^2);    
    omegac2 = omega_c*2;
    Ktemp2 = (omegac2^2./cw.^2) - repmat(kxy, 1,1,mgrid.num_z+1);  
    
    K = -sqrt(Ktemp); 
    K2 = -sqrt(Ktemp2); 
    clear Ktemp1  Ktemp2
    K(K==0) = eps;
    K2(K2==0) = eps;    
    
%     c_c = medium.c;
    if length(cw2) == 1
        cw2 = cw2*ones(size(M_linear));
    end
    
%     rho_rho = medium.rho;
%     if length(medium.rho) == 1
%         rho_rho = medium.rho*ones(size(M_linear));
%     end     
       
    % preallocate an array to store the pressure at the second harmonic
    P_second = zeros(mgrid.num_x,mgrid.num_y, mgrid.num_z+1);

    M_source = -2j.*omega_c^2.*medium.beta./medium.rho./cw.^4.*(P_fundamental).^2;
    excit_F = 0; 
    f = 0;
    
    expn2   = exp(0.5*1i*K*mgrid.dz); 
    
    h = waitbar(0,'Second-harmonic forward projection, please wait...');
    tic
    for I = 2:mgrid.num_z+1
        
        Ktemp1 = ((2*omega_c)^2./(cw2(:,:, I)).^2) - ...
                 evanescent_fiter*(mgrid.kx.'*ones(1, mgrid.num_y)).^2 - ...
                 evanescent_fiter*(ones(mgrid.num_x, 1)*mgrid.ky).^2;
 
        % Simpson  
        % 1    
        M  = fftshift(fft2(M_linear(:,:,I-1).*f)) + fftshift(fft2(M_source(:,:,I-1)));
        F1 = excit_F.*expn2 + 0.5*mgrid.dz*M.*expn2./(2i.*K);  % P01(z+mgrid.dz/2)
        F1(isnan(F1)) = 0;
        F1(real(Ktemp1)<=0) = 0; 
        % 2
        f  = ifft2(ifftshift(F1));  % p0(z+mgrid.dz/2)
        clear F1
        M1 = fftshift(fft2(M_linear(:,:,I-1).*f)) + fftshift(fft2(M_source(:,:,I-1)));
        F2 = excit_F.*expn2 + 0.25*mgrid.dz*expn2./(2i.*K).*(M + M1./expn2);% P12(z+mgrid.dz/2)
        clear M1
        F2(isnan(F2)) = 0;
        F2(real(Ktemp1)<=0) = 0;     
        % 3   
        f  = ifft2(ifftshift(F2)); % p1(z+mgrid.dz/2)  
        M2 = fftshift(fft2(M_linear(:,:,I-1).*f))+ fftshift(fft2(M_source(:,:,I-1)));
        F3 = F2.*expn2 + 0.5*mgrid.dz*M2.*expn2./(2i.*K);  % P03(z+mgrid.dz)
        F3(isnan(F3)) = 0;
        F3(real(Ktemp1)<=0) = 0;     
        % 4     
        f  = ifft2(ifftshift(F3)); % p0(z+mgrid.dz)  
        clear F3
        M3 = fftshift(fft2(M_linear(:,:,I).*f))+ fftshift(fft2(M_source(:,:,I)));
        F4 = F2.*expn2 + 0.25*mgrid.dz.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z+mgrid.dz)
        clear M3  F2
        F4(isnan(F4)) = 0;
        F4(real(Ktemp1)<=0) = 0;     
        % 5   
        f  = ifft2(ifftshift(F4)); % p0(z+mgrid.dz)  
        M4 = fftshift(fft2(M_linear(:,:,I).*f))+ fftshift(fft2(M_source(:,:,I)));     
        % P25(z+mgrid.dz)
        F5 = excit_F.*exp(1i*K*mgrid.dz) +...
             mgrid.dz/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*mgrid.dz)).*exp(1i*K*mgrid.dz)./(2i.*K); 
        clear M  M2  M4  F4
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
    
        f = ifft2(ifftshift(F5));
        excit_F = F5;
        clear F5
    
        %% add reflection 
%         rho2 = rho_rho(:,:,I);
%         if reflection_order~=0
%             c1 = c_c(:,:,I+1);
%             c2 = c_c(:,:,I);
%             rho1 = rho_rho(:,:,I+1);
%             T_rhoc = 2.*rho1.*c1./(rho1.*c1 + rho2.*c2);  %%layered medium
%             P_in = excit_F.*exp(1i.*K2(:,:,I+1)*mgrid.dz);
%             p_ref(:,:,I) = (ifft2(ifftshift(P_in))).*(T_rhoc-1); 
%             clear P_in T_rhoc  c1 c2 rho1
%         end        
        
        P_second(:,:,I) = f;
        waitbar(I/mgrid.num_z)
    
    end
    close(h)
    toc

elseif MDMC == 1
    
    % calculate the M term on RHS of solution equation
    [M_linear, Ktemp, cw2, cw] = Mterm3D_sec(mgrid, medium, omega_c);
    
    % add  non-reflecting layer
    M_linear = M_linear + repmat(M_gamma,1,1,mgrid.num_z+1);
    
    omegac2 = 2*omega_c;
    % read paper evaluation of a wave-vector...when omega is positive, K is negative. 
    kxy     = squeeze((mgrid.kx.'*ones(1, mgrid.num_y)).^2 +...
             (ones(mgrid.num_x, 1)*mgrid.ky).^2);
        
    Ktemp2 = (omegac2^2./cw.^2) - repmat(kxy, 1,1,mgrid.num_z+1);  
    
    K = -sqrt(Ktemp); 
    K2 = -sqrt(Ktemp2); 
    K(K==0) = eps;
    K2(K2==0) = eps;
    
%     c_c = medium.c;
%     if length(medium.c) == 1
%         c_c = medium.c*ones(size(M_linear));
%     end
    
    rho_rho = medium.rho;
    if length(medium.rho) == 1
        rho_rho = medium.rho*ones(size(M_linear));
    end     
    
    
    % preallocate an array to store the pressure at the second harmonic
    P_second = zeros(mgrid.num_x,mgrid.num_y, mgrid.num_z+1);

    M_source = -2j.*omega_c^2.*medium.beta./medium.rho./cw.^4.*(P_fundamental).^2;
    excit_F = 0; 
    f = 0;    
    
    expn2   = exp(0.5*1i*K*mgrid.dz);  
    
    h = waitbar(0,'Second-harmonic forward projection, please wait...');
    tic
    for I = 2:mgrid.num_z+1
        
        Ktemp1 = ((2*omega_c)^2./(cw2(:,:, I)).^2) - ...
                 evanescent_fiter*(mgrid.kx.'*ones(1, mgrid.num_y)).^2 - ...
                 evanescent_fiter*(ones(mgrid.num_x, 1)*mgrid.ky).^2;
 
        % Simpson  
        % 1    
        M  = fftshift(fft2(M_linear(:,:,I-1).*f)) + fftshift(fft2(M_source(:,:,I-1)));
        F1 = excit_F.*expn2 + 0.5*mgrid.dz*M.*expn2./(2i.*K);  % P01(z+mgrid.dz/2)
        F1(isnan(F1)) = 0;
        F1(real(Ktemp1)<=0) = 0; 
        % 2
        f  = ifft2(ifftshift(F1));  % p0(z+mgrid.dz/2)
        clear F1
        M1 = fftshift(fft2(M_linear(:,:,I).*f)) + fftshift(fft2(M_source(:,:,I-1)));
        F2 = excit_F.*expn2 + 0.25*mgrid.dz*expn2./(2i.*K).*(M + M1./expn2);% P12(z+mgrid.dz/2)
        clear M1
        F2(isnan(F2)) = 0;
        F2(real(Ktemp1)<=0) = 0;     
        % 3   
        f  = ifft2(ifftshift(F2)); % p1(z+mgrid.dz/2)  
        M2 = fftshift(fft2(M_linear(:,:,I-1).*f))+ fftshift(fft2(M_source(:,:,I-1)));
        F3 = F2.*expn2 + 0.5*mgrid.dz*M2.*expn2./(2i.*K);  % P03(z+mgrid.dz)
        F3(isnan(F3)) = 0;
        F3(real(Ktemp1)<=0) = 0;     
        % 4     
        f  = ifft2(ifftshift(F3)); % p0(z+mgrid.dz)  
        clear F3
        M3 = fftshift(fft2(M_linear(:,:,I).*f))+ fftshift(fft2(M_source(:,:,I)));
        F4 = F2.*expn2 + 0.25*mgrid.dz.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z+mgrid.dz)
        clear M3  F2
        F4(isnan(F4)) = 0;
        F4(real(Ktemp1)<=0) = 0;     
        % 5   
        f  = ifft2(ifftshift(F4)); % p0(z+mgrid.dz)  
        M4 = fftshift(fft2(M_linear(:,:,I).*f))+ fftshift(fft2(M_source(:,:,I)));     
        % P25(z+mgrid.dz)
        F5 = excit_F.*exp(1i*K*mgrid.dz) + ...
             mgrid.dz/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*mgrid.dz)).*exp(1i*K*mgrid.dz)./(2i.*K); 
        clear M  M2  M4  F4
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
    
        
%% corrections  

        c1 = cw(:,:,I-1);
        c2 = cw(:,:,I);
        rho1 = rho_rho(:,:,I-1);
        rho2 = rho_rho(:,:,I);   
               
        wx2 = omegac2.*omegac2;
        e1 = exp(1i*(K - K2(:,:,I))*mgrid.dz);
        k_dif2 = (wx2./medium.c0./medium.c0 - wx2./c2./c2)*mgrid.dz;      
        k_dif  = (wx2./medium.c0./medium.c0 - wx2./c1./c1)*mgrid.dz; 
        clear  wx2  
        correction = excit_F.*exp(1i*K2(:,:,I)*mgrid.dz).*...
                    (1-e1.*(1+k_dif./4i./K)./(1-k_dif2./4i./K));        
        F5 = F5 + correction;  
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
         
        T_rhoc = 2.*rho2.*c2./(rho1.*c1 + rho2.*c2);  % normal incidence
        f = (ifft2(ifftshift(F5)));
        f = f./T_rhoc;                
                
        F5 = fftshift(fft2(f));
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
   
        %%        
        f = ifft2(ifftshift(F5));
        excit_F = F5;
        clear F5
    
        %% add reflection 
%         if reflection_order~=0
%             P_in = excit_F.*exp(1i.*K2(:,:,I+1)*mgrid.dz);
%             p_ref(:,:, I) = (ifft2(ifftshift(P_in))).*(T_rhoc-1); 
%             clear P_in
%         end        
        
        P_second(:,:,I) = f;
        waitbar(I/mgrid.num_z)
    
    end
    close(h)
    toc
   
%%   reflection for the 2nd harmonic frequency, will be released in the future
% even reflection option is added, if medium has constant density and speed
% of sound, there is no reflection
% if (max(max(max(medium.c)))   ~= min(min(min(medium.c)))) || ...
%     (max(max(max(medium.rho))) ~= min(min(min(medium.rho))))
%     reflection_option = 1;  % flag for reflection
% else
%     reflection_option = 0;
% end    
    
% if reflection_order ~=0 && (reflection_option == 1)
%     
%     if MDMC == 1 % reflection propagation with correction
%         reflection = MReflection3D_fund(mgrid, medium, p_ref, ...
%                      M_linear, K, K2, cw, omegac2, reflection_order);
%         P_second = P_second + reflection; 
%     else    % reflection propagation without correction
%         reflection = Reflection3D_fund(mgrid, medium, p_ref, ...
%                      M_linear, K, cw, omegac2, reflection_order);
%         P_second = P_second + reflection;        
%         
%     end
%        
% end


end







