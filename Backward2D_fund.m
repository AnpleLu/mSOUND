function [P_fundamental] = Backward2D_fund(mgrid, medium, excit_p, omegac, ...
                           reflection_order, varargin)
                       
% DESCRIPTION:
% Simpson-like 2D integration in the forward direction at the fundamental
% frequency

% USAGE:
% [P_fundamental] = Forward2D_fund(mgrid, medium, excit_p, omegac, ...
%                   varargin)

% INPUTS:
% mgrid              input structure to define the computational domain
% medium             medium properties
% excit_p            excitation signals
% omegac             center frequency 2*pi*fc
% reflection_order   the maximum order of reflection included in the simulation

% OPTIONAL INOUTS:
% 'NRL'              add non-reflecting layer at the boundary to reduce the
% spatial aliasing error
% 'correction'       add correction due to strong heterogeneities
% 'no_refl_correction' turnoff the reflection correction in backward
% projection

% OUTPUTS:
% P_fundamental      fundamental pressure dirtribution through the domain

%%

% 1.05 is an empirical number, it can be higher if necessay
evanescent_fiter = 1.05;

MDMC = 0;  % flag for correction
MDMA = 0;  % flag for non-reflecting layer
MDMS = 1;
if ~isempty(varargin)
    for input_index = 1:length(varargin)
        switch varargin{input_index}
            case 'correction'
            MDMC = 1;   % flag for transmission corrections
            case 'NRL'
            MDMA = 1;   % flag for non-reflecting layer 
            case 'no_refl_correction'
            MDMS = 0;   % flag for reflection correction
        end
    end
end

%% non-reflecting layer
if MDMA ==1
    M_gamma = -medium.NRL_gamma./(cosh(medium.NRL_alpha.*mgrid.abx_vec.')).^2.*1i.*omegac;
else 
    M_gamma = 0;
end

if reflection_order~=0
    % preallocate an array to store the reflected part
    p_ref = zeros(mgrid.num_x, mgrid.num_y+1);  
end

%%
backward_dy = -mgrid.dy;
if MDMC == 0 
    % calculate the M term on RHS of solution equation
    [M_linear, Ktemp, Ktemp2, cw] = Mterm2D_fund(mgrid, medium, omegac);
    
    % add non-reflecting layer
    M_linear = M_linear + repmat(M_gamma,1,mgrid.num_y+1);
 
    if length(cw)==1
        cw = cw*ones(size(M_linear));
    end
    
    c_c = medium.c;
    if length(medium.c) == 1
        c_c = medium.c*ones(size(M_linear));
    end
    
    rho_rho = medium.rho;
    if length(medium.rho) == 1
        rho_rho = medium.rho*ones(size(M_linear));
    end    
    
    % read paper evaluation of a wave-vector...when omega is positive, K is negative.
    K  = -sqrt(Ktemp); 
    K2 = -sqrt(Ktemp2);  
    clear Ktemp Ktemp2
    K(K==0)   = eps;
    K2(K2==0) = eps; 
    
    
    % exponential term for the forward propagation
    expn2 = exp(0.5*1i*K*backward_dy); 

    % preallocate an array to store the pressure at the fundamental frequency
    P_fundamental = zeros(mgrid.num_x, mgrid.num_y+1);
    P_fundamental(:,end) = excit_p; 
     
    % normalize the pressure field
    f = excit_p.'./sqrt(rho_rho(:,end));
    
    % Fourier transform of the excitation with regard of x
    excit_F = fftshift(fft(f));     
    
    h = waitbar(0,'Backward projection, please wait...');

    for I = mgrid.num_y:-1:1
        
        Ktemp1 = (omegac).^2./cw(:,I).^2 - evanescent_fiter*mgrid.kx.'.^2;
        
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
        excit_F = F5;
        clear F5
        
        rho2 = rho_rho(:,I);
        if reflection_order~=0
            c1 = c_c(:,I-1);
            c2 = c_c(:,I);
            rho1 = rho_rho(:,I-1);            
            T_rhoc = 2.*rho1.*c1./(rho1.*c1 + rho2.*c2);  %%layered medium
            P_in = excit_F.*exp(1i.*K2(:,I+1)*mgrid.dy);
            p_ref(:,I) = (ifft(ifftshift(P_in))).*(T_rhoc-1);  
            clear T_rhoc  P_in c1  c2 rho1
        end
        % recover pressure from the normalized wave field
        P_fundamental(:,I) = f.*sqrt(rho2);
       
        waitbar(I/mgrid.num_y)
    end
    % close waitbar
    close(h)    
    
    
elseif MDMC ==1
    
    % calculate the M term on RHS of solution equation
    [M_linear, Ktemp, Ktemp2, cw] = Mterm2D_Mfund(mgrid, medium, omegac);
    
    % add non-reflecting layer
    M_linear = M_linear + M_gamma;
    
    if length(cw)==1
        cw = cw*ones(size(M_linear));
    end
    
    c_c = medium.c;
    if length(medium.c) == 1
        c_c = medium.c*ones(size(M_linear));
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
    expn2 = exp(0.5*1i*K*backward_dy); 

    % preallocate an array to store the pressure at the fundamental frequency
    P_fundamental = zeros(mgrid.num_x, mgrid.num_y+1);
    P_fundamental(:,end) = excit_p; 

    f = excit_p.';
    
    % Fourier transform of the excitation with regard of x,y and time
    excit_F = fftshift(fft(f));     
    
    h = waitbar(0,'Backward projection, please wait...');
    
    for I = mgrid.num_y:-1:1
        

        Ktemp1 = (omegac).^2./cw(:,I).^2 - evanescent_fiter*mgrid.kx.'.^2;
        
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
        %% phase correction 
        
        c1 = c_c(:,I+1);
        c2 = c_c(:,I);
        rho1 = rho_rho(:,I+1);
        rho2 = rho_rho(:,I);

        wx2 = omegac.*omegac;
        e1 = exp(1i*(K - K2(:,I))*backward_dy);
        k_dif2 = (wx2./medium.c0./medium.c0 - wx2./c2./c2)*backward_dy;      
        k_dif  = (wx2./medium.c0./medium.c0 - wx2./c1./c1)*backward_dy; 
        clear  wx2  
        correction = excit_F.*exp(1i*K2(:,I)*backward_dy).*...
                    (1-e1.*(1+k_dif./4i./K)./(1-k_dif2./4i./K));
        F5 = F5 + correction;  
        clear correction   el  k_dif2  k_dif
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
  %%   transmission
        T_rhoc = 2.*rho1.*c1./(rho1.*c1 + rho2.*c2); 
        f = (ifft(ifftshift(F5)));
        f = f./T_rhoc;       
 %%   
        F5 = fftshift(fft(f));
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
        f = ifft(ifftshift(F5));
      
        excit_F = F5;
        clear F5
    
        if reflection_order~=0
            P_in = excit_F.*exp(1i.*K2(:,I+1)*mgrid.dy);
            p_ref(:,I) = (ifft(ifftshift(P_in))).*(T_rhoc-1);  
            clear  P_in
        end        
       P_fundamental(:,I) = f;
       waitbar(I/mgrid.num_y)
    end
    
    close(h)

end

% even reflection option is added, it medium has constant density and speed
% of sound, there is no reflection
if (max(max(max(medium.c)))   ~= min(min(min(medium.c)))) || ...
    (max(max(max(medium.rho))) ~= min(min(min(medium.rho))))
    reflection_option = 1;  % flag for reflection
else
    reflection_option = 0;
end

if reflection_order ~=0 && (reflection_option == 1)
    
    if MDMC == 1 % reflection propagation with correction
        reflection = BMReflection2D_fund(mgrid, medium, p_ref, ...
                     M_linear, K, K2, cw, omegac,...
                     reflection_order, excit_p, MDMS);
                 
        P_fundamental = P_fundamental + reflection; 
    else    % reflection propagation without correction
        reflection = BReflection2D_fund(mgrid, medium, p_ref, ...
                     M_linear, K, cw, omegac,...
                     reflection_order, excit_p, MDMS);
                 
        P_fundamental = P_fundamental + reflection;        
        
    end

    if MDMS == 1  % with reflection correction
        P_fundamental = reflection; 
    else    % without reflection correction
        P_fundamental = P_fundamental + reflection;
    end    
     
       
end

end