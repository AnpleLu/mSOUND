function [P_fundamental] = Forward2D_fund(mgrid, medium, excit_p, omega_c, ...
                           reflection_order, varargin)
                       
% DESCRIPTION:
% Simpson-like 2D integration in the forward direction at the fundamental
% frequency

% USAGE:
% [P_fundamental] = Forward2D_fund(mgrid, medium, excit_p, omega_c, ...
%                   reflection_order, varargin)

% INPUTS:
% mgrid              input structure to define the computational domain
% medium             medium properties
% excit_p            excitation signals
% omega_c            center frequency 2*pi*fc
% reflection_order   the maximum order of reflection included in the simulation

% OPTIONAL INOUTS:
% 'NRL'              add non-reflecting layer at the boundary to reduce the
% spatial aliasing error
% 'correction'       add correction due to strong heterogeneities

% OUTPUTS:
% P_fundamental      fundamental pressure dirtribution through the domain

%%
MDMC = 0;  % flag for correction
MDMA = 0;  % flag for non-reflecting layer
if ~isempty(varargin)
    for input_index = 1:length(varargin)
        switch varargin{input_index}
            case 'correction'
            MDMC = 1;   % flag for transmission correction
            case 'NRL'
            MDMA = 1;   % flag for non-reflection layer
        end
    end
end

%% non-reflecting layer
if MDMA ==1
    M_gamma = medium.NRL_gamma./(cosh(medium.NRL_alpha.*mgrid.abx_vec.')).^2.*1i.*omega_c;
else 
    M_gamma = 0;
end 

% 1.05 is an empirical number, it can be higher if necessay
evanescent_fiter = 1.05;

% figure
% plot(mgrid.x*1e3, abs(M_gamma), 'r','LineWidth',2)
%%
if reflection_order~=0
    % preallocate an array to store the reflected part
    p_ref = zeros(mgrid.num_x, mgrid.num_y+1);  
end
    
if MDMC == 0 
    % calculate the M term on RHS of solution equation
    [M_linear, Ktemp, Ktemp2, cw] = Mterm2D_fund(mgrid, medium, omega_c);
    
    % add non-reflecting layer
    M_linear = M_linear + repmat(M_gamma,1,mgrid.num_y+1);
    
    if length(cw)==1
        cw = cw*ones(size(M_linear));
    end

    rho_rho = medium.rho;
    if length(medium.rho) == 1
        rho_rho = medium.rho*ones(size(M_linear));
    end
    % read paper evaluation of a wave-vector...when omega is positive, K is negative.
    K  = -sqrt(Ktemp); 
    K2 = -sqrt(Ktemp2); 
    K(K==0)   = eps;
    K2(K2==0) = eps;
    
    % exponential term for the forward propagation
    expn2 = exp(0.5*1i*K*mgrid.dy); 

    % preallocate an array to store the pressure at the fundamental frequency
    P_fundamental = zeros(mgrid.num_x, mgrid.num_y+1);
    P_fundamental(:,1) = excit_p;
    
    % normalize the pressure field
    f = excit_p.'./sqrt(rho_rho(:,1));
    
    % Fourier transform of the excitation with regard of x
    excit_F = fftshift(fft(f));     
    
    h = waitbar(0,'Forward projection, please wait...');

    for I = 2:mgrid.num_y+1
        
        Ktemp1 = (omega_c).^2./cw(:,I).^2 -evanescent_fiter*mgrid.kx.'.^2;
        
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
        F5 = excit_F.*exp(1i*K*mgrid.dy) +...
             mgrid.dy/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*mgrid.dy)).*exp(1i*K*mgrid.dy)./(2i.*K); 
        clear M  M2  M4  F4
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
        
        %% add reflection 
        c1 = cw(:,I-1);
        c2 = cw(:,I);
        rho1 = rho_rho(:,I-1);
        rho2 = rho_rho(:,I);
        
        if reflection_order~=0
            T_rhoc = 2.*rho2.*c2./(rho1.*c1 + rho2.*c2);  %%layered medium
            P_in = excit_F.*exp(1i.*K2(:,I)*mgrid.dy);
            p_ref(:,I-1) = (ifft(ifftshift(P_in))).*(T_rhoc-1);  
            clear T_rhoc  P_in
        end
          
        excit_F = F5;
    
        f  = ifft(ifftshift(F5)); % p0(z+mgrid.dy) 
         clear F5
        % recover pressure from the normalized wave field
        P_fundamental(:,I) = f.*sqrt(rho2);
       
        waitbar(I/mgrid.num_y)
    end
    % close waitbar
    close(h)    
    
    
elseif MDMC ==1
    
    % calculate the M term on RHS of solution equation
    [M_linear, Ktemp, Ktemp2, cw] = Mterm2D_Mfund(mgrid, medium, omega_c);
    
    % add non-reflecting layer
    M_linear = M_linear + M_gamma;
    
    if length(cw)==1
        cw = cw*ones(size(M_linear));
    end
    
    rho_rho = medium.rho;
    if length(medium.rho) == 1
        rho_rho = medium.rho*ones(size(M_linear));
    end   
    % read paper evaluation of a wave-vector...when omega is positive, K is negative.
    K = -sqrt(Ktemp); 
    K2 = -sqrt(Ktemp2); 
    K(K==0)   = eps;
    K2(K2==0) = eps;
    clear  Ktemp   Ktemp2
    
    % exponential term for the forward propagation
    expn2 = exp(0.5*1i*K*mgrid.dy); 

    % preallocate an array to store the pressure at the fundamental frequency
    P_fundamental = zeros(mgrid.num_x, mgrid.num_y+1);
    P_fundamental(:,1) = excit_p;

    f = excit_p.';
    
    % Fourier transform of the excitation with regard of x,y and time
    excit_F = fftshift(fft(f));     
    
    h = waitbar(0,'Forward projection, please wait...');
    
    for I = 2:mgrid.num_y+1
        
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
        F5 = excit_F.*exp(1i*K*mgrid.dy) +...
             mgrid.dy/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*mgrid.dy)).*exp(1i*K*mgrid.dy)./(2i.*K); 
        clear M  M2  M4  F4
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
        %% phase correction 
        c1 = cw(:,I-1);
        c2 = cw(:,I);
        rho1 = rho_rho(:,I-1);
        rho2 = rho_rho(:,I);

        wx2 = omega_c.*omega_c;
        e1 = exp(1i*(K - K2(:,I))*mgrid.dy);
        k_dif2 = (wx2./medium.c0./medium.c0 - wx2./c2./c2)*mgrid.dy;      
        k_dif  = (wx2./medium.c0./medium.c0 - wx2./c1./c1)*mgrid.dy; 
        clear  wx2  
        correction = excit_F.*exp(1i*K2(:,I)*mgrid.dy).*...
                    (1-e1.*(1+k_dif./4i./K)./(1-k_dif2./4i./K));
        F5 = F5 + correction;  
        clear correction   el  k_dif2  k_dif
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
  %%   transmission
        T_rhoc = 2.*rho2.*c2./(rho1.*c1 + rho2.*c2); 
        f = (ifft(ifftshift(F5)));
        f = f.*T_rhoc;       
 %%   
        F5 = fftshift(fft(f));
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
        %% calculate reflected part
        if reflection_order~=0
            P_in = excit_F.*exp(1i.*K2(:,I)*mgrid.dy);
            p_ref(:,I-1) = (ifft(ifftshift(P_in))).*(T_rhoc-1);  
            clear P_in
        end
        
        f = ifft(ifftshift(F5));
        excit_F = F5;
        clear F5
    
        P_fundamental(:,I) = f;
        waitbar(I/mgrid.num_y)
    end
    
    close(h)

end

% even reflection option is added, it medium has constant density and speed
% of sound, there is no reflection
if (max(max(medium.c))   ~= min(min(medium.c))) || ...
    (max(max(medium.rho)) ~= min(min(medium.rho)))
    reflection_option = 1;  % flag for reflection
else
    reflection_option = 0;
end
    
if reflection_order ~=0 && (reflection_option == 1)
    
    if MDMC == 1 % reflection propagation with correction
        reflection = MReflection2D_fund(mgrid, medium, p_ref, ...
                     M_linear, K, K2, cw, omega_c, reflection_order);
        P_fundamental = P_fundamental + reflection; 
    else    % reflection propagation without correction
        reflection = Reflection2D_fund(mgrid, medium, p_ref, ...
                     M_linear, K, cw, omega_c, reflection_order);
        P_fundamental = P_fundamental + reflection;        
        
    end
       
end


end