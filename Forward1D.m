function [p_time] = Forward1D(mgrid, medium, excit_ps, sensor_mask,...
                    reflection_order, varargin)
% DESCRIPTION:
% Simpson-like integration in the forward direction 

% USAGE:
% [p_time] = Forward1D(mgrid, medium, excit_ps, sensor_mask,...
%                      reflection_order, varargin)

% INPUTS:
% mgrid              structure to define the computational domain
% medium             medium properties
% excit_ps           excitation signals
% sensor_mask        a set of Cartesian points where the pressure is recorded
% reflection_order   the maximum order of reflection included in the simulation

% OPTIONAL INOUTS
% 'correction'       apply phase and amplitude corrections for strongly
%                    heterogeneous media

% OUTPUTS:
% p_time             total pressure dirtribution through the domain

%%
% reshape the excitation signal for forward projection
excit_p = excitation(excit_ps, mgrid);
clear excit_ps

% define p_time to store the time-domain results
p_time = zeros(mgrid.num_t, sum(sum(sensor_mask)));

if sum(sensor_mask(1))~=0
    p_time(:,1) = excit_p;
end

% define p_ref to store the reflected part
p_ref = zeros(mgrid.num_t,mgrid.num_x+1);

% find medium properties
medium_cond = medium_case(medium);
MDMC = 0;

if ~isempty(varargin)
    for input_index = 1:length(varargin)
        switch varargin{input_index}
            case 'correction'
            MDMC = 1;     % flag for correction
        end
    end
end


%%
if medium_cond == 1
    
    disp( 'linear homogeneous medium')
    
     % calculate the M term on RHS of solution equation
    [~, M_linear, cw] = Mterm1D(mgrid, medium, medium_cond);
     

     Ktemp = mgrid.w.'.^2./cw.^2;  % with constant c   
     K = sqrt(Ktemp);
     
     if mod (mgrid.num_t,2) == 1
         K(mgrid.num_t/2+0.5:end,:)  = -K(mgrid.num_t/2+0.5:end,:);
     else
         K(mgrid.num_t/2:end,:)  = -K(mgrid.num_t/2:end,:);
     end
     K(K==0) = eps;   
        
     f = excit_p;
     clear  excit_p  Ktemp
     
     % exponential term for the forward propagation
     expn2 = exp(0.5*1i*K*mgrid.dx); 
     
     % Fourier transform of the excitation with regard of time
     excit_F = fftshift(fft(f)); 
     
     h = waitbar(0,'Forward projection, please wait...'); 
     
     for I = 2:mgrid.num_x+1  
         
         %Integral algorithm
         % 1   
         Ft  = fftshift(fft(f,[],1),1);
         % M(f(z))
         M = M_linear.*Ft;
         F1 = excit_F.*expn2 + 0.5*mgrid.dx*M.*expn2./(2i.*K);  % P01(z+mgrid.dx/2)
         F1(isnan(F1)) = 0;
    
         % 2
         f = real(ifft(ifftshift(F1)));  % p0(z+mgrid.dx/2) 
         clear F1
         Ft  = fftshift(fft(f,[],1),1);     % fft(f)   
         % M(f0(z+mgrid.dx/2))  M1
         M1 = M_linear.*Ft;
         F2 = excit_F.*expn2 + 0.25*mgrid.dx*expn2./(2i.*K).*(M + M1./expn2);% P12(z+mgrid.dx/2)
         F2(isnan(F2)) = 0;
   
         % 3   
         f   = real(ifft(ifftshift(F2))); % p1(z+mgrid.dx/2)   
         Ft  = fftshift(fft(f,[],1),1);    % fft(f1) 
         % M(f1(z+mgrid.dx/2)) M2
         M2 = M_linear.*Ft;
         F3 = F2.*expn2 + 0.5*mgrid.dx*M2.*expn2./(2i.*K);  % P03(z+mgrid.dx)
         F3(isnan(F3)) = 0;
    
         % 4     
         f   = real(ifft(ifftshift(F3))); % p0(z+mgrid.dx)  
         clear F3
         Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
         % M(f0(z+mgrid.dx))  M3
         M3 = M_linear.*Ft;
         F4 = F2.*expn2 + 0.25*mgrid.dx.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z+mgrid.dx)
         clear F2
         F4(isnan(F4)) = 0;  
      
         % 5   
         f   = real(ifft(ifftshift(F4))); % p0(z+mgrid.dx)
         clear F4
         Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
         % M(f1(z+mgrid.dx)) M4
         M4 = M_linear.*Ft;
         F5 = excit_F.*exp(1i*K*mgrid.dx) +...
              mgrid.dx/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*mgrid.dx)).*...
              exp(1i*K*mgrid.dx)./(2i.*K); % P25(z+mgrid.dx)
     
         clear M  M2  M4  
         F5(isnan(F5)) = 0;

         excit_F = F5;
         f = real(ifft(ifftshift(F5))); % p0(z+mgrid.dx) 
         clear F5 
    
         if sum(sensor_mask(I))~=0
            
             % calculate the starting position of recorded signal
             sensor_mask_I1 = sum((sensor_mask(1:I-1)))+1;
        
             % calculate the ending position of recorded signal
             sensor_mask_I2 = sensor_mask_I1 -1 + ...
             sum(sensor_mask(I));
        
             % save the time-domain signals
             p_time(:,sensor_mask_I1:sensor_mask_I2) = f;
         end
     
         waitbar(I/mgrid.num_x)
     end
     close(h)    
 
elseif medium_cond == 2 || medium_cond == 3
    
    disp ('Inhomogeneous only in nonlinearity coefficient')
        
    % normalized wave field
    f = excit_p;
    clear excit_p    
    
    % Fourier transform of the excitation with regard of x and time
    excit_F = fftshift(fft(f)); 
        
    % calculate the M term on RHS of solution equation
    [M_nonlinear, M_linear, cw] = Mterm1D(mgrid, medium, medium_cond);

     Ktemp = mgrid.w.'.^2./cw.^2;  % with constant c   
     K = sqrt(Ktemp);
   
     if mod (mgrid.num_t,2) == 1
         K(mgrid.num_t/2+0.5:end,:)  = -K(mgrid.num_t/2+0.5:end,:);
     else
         K(mgrid.num_t/2:end,:)  = -K(mgrid.num_t/2:end,:);
     end
     K(K==0) = eps;  
     clear Ktemp
        
     % exponential term for the forward propagation
     expn2 = exp(0.5*1i*K*mgrid.dx); 
    
     % waitbar
     h = waitbar(0,'Forward projection, please wait...');

     for I = 2:mgrid.num_x+1  
         
         %Simpson integral algorithm

         % 1   
         Ft  = fftshift(fft(f,[],1),1);
         F2t = fftshift(fft(f.*f,[],1),1);
         % M(f(z))
         M = M_linear.*Ft + M_nonlinear(:,I-1).*F2t;
         F1 = excit_F.*expn2 + 0.5*mgrid.dx*M.*expn2./(2i.*K);  % P01(z+mgrid.dx/2)
         F1(isnan(F1)) = 0;
    
         % 2
         f = real(ifft(ifftshift(F1)));  % p0(z+mgrid.dx/2) 
         clear F1
         Ft  = fftshift(fft(f,[],1),1);     % fft(f) 
         F2t = fftshift(fft(f.*f,[],1),1);  % fft(f^2)   
         % M(f0(z+mgrid.dx/2))  M1
         M1 = M_linear.*Ft + M_nonlinear(:,I-1).*F2t;
         F2 = excit_F.*expn2 + 0.25*mgrid.dx*expn2./(2i.*K).*(M + M1./expn2);% P12(z+mgrid.dx/2)
         F2(isnan(F2)) = 0;
   
         % 3   
         f   = real(ifft(ifftshift(F2))); % p1(z+mgrid.dx/2)   
         Ft  = fftshift(fft(f,[],1),1);    % fft(f1) 
         F2t = fftshift(fft(f.*f,[],1),1); % fft(f1^2)   
         % M(f1(z+mgrid.dx/2)) M2
         M2 = M_linear.*Ft + M_nonlinear(:,I-1).*F2t;
         F3 = F2.*expn2 + 0.5*mgrid.dx*M2.*expn2./(2i.*K);  % P03(z+mgrid.dx)
         F3(isnan(F3)) = 0;
    
         % 4     
         f   = real(ifft(ifftshift(F3))); % p0(z+mgrid.dx)  
         clear F3
         Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
         F2t = fftshift(fft(f.*f,[],1),1); % fft(f2^2)
         % M(f0(z+mgrid.dx))  M3
         M3 = M_linear.*Ft + M_nonlinear(:,I).*F2t;
         F4 = F2.*expn2 + 0.25*mgrid.dx.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z+mgrid.dx)
         clear F2
         F4(isnan(F4)) = 0;  
      
         % 5   
         f   = real(ifft(ifftshift(F4))); % p0(z+mgrid.dx)
         clear F4
         Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
         F2t = fftshift(fft(f.*f,[],1),1); % fft(f2^2)
         % M(f1(z+mgrid.dx)) M4
         M4 = M_linear.*Ft + M_nonlinear(:,I).*F2t;
         F5 = excit_F.*exp(1i*K*mgrid.dx) +...
              mgrid.dx/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*mgrid.dx)).*...
              exp(1i*K*mgrid.dx)./(2i.*K); % P25(z+mgrid.dx)
     
         clear M  M2  M4  
         F5(isnan(F5)) = 0;

         excit_F = F5;
         f = real(ifft(ifftshift(F5))); % p0(z+mgrid.dx) 
         clear F5 
    
         if sum(sensor_mask(I))~=0
            
             % calculate the starting position of recorded signal
             sensor_mask_I1 = sum((sensor_mask(1:I-1)))+1;
        
             % calculate the ending position of recorded signal
             sensor_mask_I2 = sensor_mask_I1 -1 + ...
             sum(sensor_mask(I));
        
             % save the time-domain signals
             p_time(:,sensor_mask_I1:sensor_mask_I2) = f;
         end
     
         waitbar(I/mgrid.num_x)
     end
     close(h)    
     
elseif medium_cond == 4 && MDMC == 0
    
    disp ('Linear medium without transmission corrections')
    
    % calculate the M term on RHS of solution equation
    [~, M_linear, cw] = Mterm1D(mgrid, medium, medium_cond);
    
   
    % wave number along the propagation direction and it is w/c for 1D case  
    Ktemp = mgrid.w.'.^2./medium.c0.^2;    % with constant c 
    Ktemp2 = mgrid.w.'.^2./cw.^2;          % with heterogeneous c

    K  = sqrt(Ktemp);
    K2 = sqrt(Ktemp2);
    if mod (mgrid.num_t,2) == 1
        K(mgrid.num_t/2+0.5:end,:)  = -K(mgrid.num_t/2+0.5:end,:);
        K2(mgrid.num_t/2+0.5:end,:) = -K2(mgrid.num_t/2+0.5:end,:);
    else
        K(mgrid.num_t/2:end,:)  = -K(mgrid.num_t/2:end,:);
        K2(mgrid.num_t/2:end,:) = -K2(mgrid.num_t/2:end,:);
    end         
    K(K==0) = eps;
    K2(K2==0) = eps;    
    
    % exponential term for the forward propagation
    expn2 = exp(0.5*1i*K*mgrid.dx); 

    rho_rho = medium.rho;
    if length(medium.rho) == 1
        rho_rho = medium.rho*ones(mgrid.num_x+1);
    end     
    
    % normalized wave field
    f = excit_p/sqrt(rho_rho(1));
        
    % Fourier transform of the excitation with regard of x and time
    excit_F = fftshift(fft(f)); 

    % waitbar
    h = waitbar(0,'Forward projection, please wait...');

    for I = 2:mgrid.num_x+1  
        
        % Simpson integral algorithm
        c1 = cw(:,I-1);
        c2 = cw(:,I);    
        rho1 = rho_rho(I-1);
        rho2 = rho_rho(I);
        
        %Simpson integral algorithm

        % 1   
        Ft  = fftshift(fft(f,[],1),1);
        % M(f(z))
        M = M_linear(:,I-1).*Ft;
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
   
        %% reflection
        if reflection_order~=0
            T_rhoc = 2.*rho2.*c2./(rho1.*c1 + rho2.*c2);  %%layered medium   
            P_in = excit_F.*exp(1i.*K2(:,I-1)*mgrid.dx);
            p_ref(:,I-1) = real(ifft2(ifftshift(P_in))).*(T_rhoc-1);   
        end
        
        excit_F = F5;
        f = real(ifft(ifftshift(F5))); % p0(z+mgrid.dx) 
        clear F5

        % recover the normalized wave field
        p_I = sqrt(medium.rho(I)).*f;    
    
        if sum(sensor_mask(I))~=0
            
            % calculate the starting position of recorded signal
            sensor_mask_I1 = sum((sensor_mask(1:I-1)))+1;
        
            % calculate the ending position of recorded signal
            sensor_mask_I2 = sensor_mask_I1 -1 + ...
            sum(sensor_mask(I));
        
            % save the time-domain signals
            p_time(:,sensor_mask_I1:sensor_mask_I2) = p_I;
        end
     
        clear p_I
        waitbar(I/mgrid.num_x)
    end
    close(h)     
     
elseif medium_cond == 4 && MDMC == 1 % add correction
    
    disp ('Linear medium with transmission corrections')
        
    % calculate the M term on RHS of solution equation
    [~, M_linear, cw] = Mterm1D_MMDM(mgrid, medium, medium_cond);    

    % wave number along the propagation direction and it is w/c for 1D case  
    Ktemp = mgrid.w.'.^2./medium.c0.^2;    % with constant c 
    Ktemp2 = mgrid.w.'.^2./cw.^2;          % with heterogeneous c

    K  = sqrt(Ktemp);
    K2 = sqrt(Ktemp2);
    if mod (mgrid.num_t,2) == 1
        K(mgrid.num_t/2+0.5:end,:)  = -K(mgrid.num_t/2+0.5:end,:);
        K2(mgrid.num_t/2+0.5:end,:) = -K2(mgrid.num_t/2+0.5:end,:);
    else
        K(mgrid.num_t/2:end,:)  = -K(mgrid.num_t/2:end,:);
        K2(mgrid.num_t/2:end,:) = -K2(mgrid.num_t/2:end,:);
    end

    K(K==0) = eps;
    K2(K2==0) = eps;   
    
    f = excit_p;
    clear excit_p
    
    rho_rho = medium.rho;
    if length(medium.rho) == 1
        rho_rho = medium.rho*ones(mgrid.num_x+1);
    end        
    
    % Fourier transform of the excitation with regard of x and time
    excit_F = fftshift(fft(f)); 

    % exponential term for the forward propagation
    expn2 = exp(0.5*1i*K*mgrid.dx); 
    
    % waitbar
    h = waitbar(0,'Forward projection, please wait...');

    for I = 2:mgrid.num_x+1      
        % Simpson integral algorithm
        c1 = cw(:,I-1);
        c2 = cw(:,I);    
        rho1 = rho_rho(I-1);
        rho2 = rho_rho(I);
 
        % 1   
        Ft  = fftshift(fft(f,[],1),1);
        % M(f(z))
        M = M_linear(:,I-1).*Ft;
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
%% reflection        
        if reflection_order~=0
            P_in = excit_F.*exp(1i.*K2(:,I-1)*mgrid.dx);
            p_ref(:,I-1) = real(ifft2(ifftshift(P_in))).*(T_rhoc-1);   
        end   
        
        
        
        excit_F = F5;
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
            p_time(:,sensor_mask_I1:sensor_mask_I2) = f;
        
        end
        
        waitbar(I/mgrid.num_x)
    
    end
    close(h)     
     
 %%        
elseif medium_cond == 5 && MDMC == 0
    
%     disp ('Nonlinear medium without transmission corrections')
    
    % M_term 
    [M_nonlinear, M_linear, cw] = Mterm1D(mgrid, medium, medium_cond);    
    
    % wave number along the propagation direction and it is w/c for 1D case  
    Ktemp = mgrid.w.'.^2./medium.c0.^2;  % with constant c   
    Ktemp2 = mgrid.w.'.^2./cw.^2;        % with heterogeneous c

    K  = sqrt(Ktemp);
    K2 = sqrt(Ktemp2);
    
    clear Ktemp Ktemp2
    if mod (mgrid.num_t,2) == 1
        K(mgrid.num_t/2+0.5:end,:)  = -K(mgrid.num_t/2+0.5:end,:);
        K2(mgrid.num_t/2+0.5:end,:) = -K2(mgrid.num_t/2+0.5:end,:);
    else
        K(mgrid.num_t/2:end,:)  = -K(mgrid.num_t/2:end,:);
        K2(mgrid.num_t/2:end,:) = -K2(mgrid.num_t/2:end,:);
    end

    K(K==0) = eps;
    K2(K2==0) = eps;       

    rho_rho = medium.rho;
    if length(medium.rho) == 1
        rho_rho = medium.rho*ones(mgrid.num_x+1);
    end        
    
    % normalized wave field
    f = excit_p/sqrt(rho_rho(1));

    % Fourier transform of the excitation with regard of x and time
    excit_F = fftshift(fft(f));

    % exponential term for the forward propagation
    expn2 = exp(0.5*1i*K*mgrid.dx);
 
    % waitbar
    h = waitbar(0,'Forward projection, please wait...');

    for I = 2:mgrid.num_x+1  
        
        % Simpson integral algorithm
        c1 = cw(:,I-1);
        c2 = cw(:,I);    
        rho1 = rho_rho(I-1);
        rho2 = rho_rho(I);
        
        %Simpson integral algorithm

        % 1   
        Ft  = fftshift(fft(f,[],1),1);
        F2t = fftshift(fft(f.*f,[],1),1);
        % M(f(z))
        M = M_linear(:,I-1).*Ft + M_nonlinear(:,I-1).*F2t;
        F1 = excit_F.*expn2 + 0.5*mgrid.dx*M.*expn2./(2i.*K);  % P01(z+mgrid.dx/2)
        F1(isnan(F1)) = 0;
    
        % 2
        f = real(ifft(ifftshift(F1)));  % p0(z+mgrid.dx/2) 
        clear F1
        Ft  = fftshift(fft(f,[],1),1);     % fft(f) 
        F2t = fftshift(fft(f.*f,[],1),1);  % fft(f^2)   
        % M(f0(z+mgrid.dx/2))  M1
        M1 = M_linear(:,I-1).*Ft + M_nonlinear(:,I-1).*F2t;
        F2 = excit_F.*expn2 + 0.25*mgrid.dx*expn2./(2i.*K).*(M + M1./expn2);% P12(z+mgrid.dx/2)
        F2(isnan(F2)) = 0;
   
        % 3   
        f   = real(ifft(ifftshift(F2))); % p1(z+mgrid.dx/2)   
        Ft  = fftshift(fft(f,[],1),1);    % fft(f1) 
        F2t = fftshift(fft(f.*f,[],1),1); % fft(f1^2)   
        % M(f1(z+mgrid.dx/2)) M2
        M2 = M_linear(:,I-1).*Ft + M_nonlinear(:,I-1).*F2t;
        F3 = F2.*expn2 + 0.5*mgrid.dx*M2.*expn2./(2i.*K);  % P03(z+mgrid.dx)
        F3(isnan(F3)) = 0;
    
        % 4     
        f   = real(ifft(ifftshift(F3))); % p0(z+mgrid.dx)  
        clear F3
        Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
        F2t = fftshift(fft(f.*f,[],1),1); % fft(f2^2)
        % M(f0(z+mgrid.dx))  M3
        M3 = M_linear(:,I).*Ft + M_nonlinear(:,I).*F2t;
        F4 = F2.*expn2 + 0.25*mgrid.dx.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z+mgrid.dx)
        clear F2
        F4(isnan(F4)) = 0;  
      
        % 5   
        f   = real(ifft(ifftshift(F4))); % p0(z+mgrid.dx)
        clear F4
        Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
        F2t = fftshift(fft(f.*f,[],1),1); % fft(f2^2)
        % M(f1(z+mgrid.dx)) M4
        M4 = M_linear(:,I).*Ft + M_nonlinear(:,I).*F2t;
        F5 = excit_F.*exp(1i*K*mgrid.dx) +...
             mgrid.dx/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*mgrid.dx)).*...
             exp(1i*K*mgrid.dx)./(2i.*K); % P25(z+mgrid.dx)
     
        clear M  M2  M4  
        F5(isnan(F5)) = 0;
   
        %% reflection
        if reflection_order~=0
            T_rhoc = 2.*rho2.*c2./(rho1.*c1 + rho2.*c2);  %%layered medium   
            P_in = excit_F.*exp(1i.*K2(:,I-1)*mgrid.dx);
            p_ref(:,I-1) = real(ifft2(ifftshift(P_in))).*(T_rhoc-1);   
        end        
                  
        excit_F = F5;
        f = real(ifft(ifftshift(F5))); % p0(z+mgrid.dx) 
        clear F5

        % recover the normalized wave field
        p_I = sqrt(rho_rho(I)).*f;    
    
        if sum(sensor_mask(I))~=0
            
            % calculate the starting position of recorded signal
            sensor_mask_I1 = sum((sensor_mask(1:I-1)))+1;
        
            % calculate the ending position of recorded signal
            sensor_mask_I2 = sensor_mask_I1 -1 + ...
            sum(sensor_mask(I));
        
            % save the time-domain signals
            p_time(:,sensor_mask_I1:sensor_mask_I2) = p_I;
        end
     
        clear p_I
        waitbar(I/mgrid.num_x)
    end
    close(h)


%%
elseif medium_cond == 5 && MDMC == 1
    
%    disp ('Nonlinear medium with transmission corrections')
    % calculate the M term on RHS of solution equation 
    [M_nonlinear, M_linear, cw] = Mterm1D_MMDM(mgrid, medium, medium_cond);

    % wave number along the propagation direction and it is w/c for 1D case  
    Ktemp = mgrid.w.'.^2./medium.c0.^2;    % with constant c 
    Ktemp2 = mgrid.w.'.^2./cw.^2;          % with heterogeneous c

    K  = sqrt(Ktemp);
    K2 = sqrt(Ktemp2);
    if mod (mgrid.num_t,2) == 1
        K(mgrid.num_t/2+0.5:end,:)  = -K(mgrid.num_t/2+0.5:end,:);
        K2(mgrid.num_t/2+0.5:end,:) = -K2(mgrid.num_t/2+0.5:end,:);
    else
        K(mgrid.num_t/2:end,:)  = -K(mgrid.num_t/2:end,:);
        K2(mgrid.num_t/2:end,:) = -K2(mgrid.num_t/2:end,:);
    end

    K(K==0) = eps;
    K2(K2==0) = eps;    
    
    rho_rho = medium.rho;
    if length(medium.rho) == 1
        rho_rho = medium.rho*ones(mgrid.num_x+1);
    end       
    
    % exponential term for the forward propagation
    expn2 = exp(0.5*1i*K*mgrid.dx);
    
    % normalized wave field
    f = excit_p;
    clear excit_p
    
    % Fourier transform of the excitation with regard of x and time
    excit_F = fftshift(fft(f));
    
    % waitbar
    h = waitbar(0,'Forward projection, please wait...');
    for I = 2:mgrid.num_x+1      
        % Integral algorithm
        c1 = cw(:,I-1);
        c2 = cw(:,I);    
        rho1 = rho_rho(I-1);
        rho2 = rho_rho(I);
 
        % 1   
        Ft  = fftshift(fft(f,[],1),1);
        F2t = fftshift(fft(f.*f,[],1),1);
        % M(f(z))
        M = M_linear(:,I-1).*Ft + M_nonlinear(:,I-1).*F2t;
        F1 = excit_F.*expn2 + 0.5*mgrid.dx*M.*expn2./(2i.*K);  % P01(z+mgrid.dx/2)
        F1(isnan(F1)) = 0;

        % 2
        f = real(ifft(ifftshift(F1)));  % p0(z+mgrid.dx/2) 
        clear F1
        Ft  = fftshift(fft(f,[],1),1);     % fft(f) 
        F2t = fftshift(fft(f.*f,[],1),1);  % fft(f^2)   
        % M(f0(z+mgrid.dx/2))  M1
        M1 = M_linear(:,I-1).*Ft + M_nonlinear(:,I-1).*F2t;
        F2 = excit_F.*expn2 + 0.25*mgrid.dx*expn2./(2i.*K).*(M + M1./expn2);% P12(z+mgrid.dx/2)
        F2(isnan(F2)) = 0;
  
        % 3   
        f   = real(ifft(ifftshift(F2))); % p1(z+mgrid.dx/2)   
        Ft  = fftshift(fft(f,[],1),1);    % fft(f1) 
        F2t = fftshift(fft(f.*f,[],1),1); % fft(f1^2)   
        % M(f1(z+mgrid.dx/2)) M2
        M2 = M_linear(:,I-1).*Ft + M_nonlinear(:,I-1).*F2t;
        F3 = F2.*expn2 + 0.5*mgrid.dx*M2.*expn2./(2i.*K);  % P03(z+mgrid.dx)
        F3(isnan(F3)) = 0;
    
        % 4     
        f   = real(ifft(ifftshift(F3))); % p0(z+mgrid.dx)  
        clear F3
        Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
        F2t = fftshift(fft(f.*f,[],1),1); % fft(f2^2)
        % M(f0(z+mgrid.dx))  M3
        M3 = M_linear(:,I).*Ft + M_nonlinear(:,I).*F2t;
        F4 = F2.*expn2 + 0.25*mgrid.dx.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z+mgrid.dx)
        clear F2
        F4(isnan(F4)) = 0;  
      
        % 5   
        f   = real(ifft(ifftshift(F4))); % p0(z+mgrid.dx)
        clear F4
        Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
        F2t = fftshift(fft(f.*f,[],1),1); % fft(f2^2)
        % M(f1(z+mgrid.dx)) M4
        M4 = M_linear(:,I).*Ft + M_nonlinear(:,I).*F2t;
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
%% reflection 
        if reflection_order~=0
            P_in = excit_F.*exp(1i.*K2(:,I-1)*mgrid.dx);
            p_ref(:,I-1) = real(ifft2(ifftshift(P_in))).*(T_rhoc-1);   
        end           
   
        excit_F = F5;
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
            p_time(:,sensor_mask_I1:sensor_mask_I2) = f;
        
        end
        
        waitbar(I/mgrid.num_x)
    
    end
    close(h)
    
end


% flag for reflection
if (max(max(max(medium.c)))   ~= min(min(min(medium.c)))) || ...
    (max(max(max(medium.rho))) ~= min(min(min(medium.rho))))
    reflection_option = 1; 
else
    reflection_option = 0;
end
    
% calcuate the reflections
if reflection_order ~=0 && (reflection_option ==1)
    if MDMC == 1 % reflection with correction
        reflection = MReflection1D(mgrid, medium, p_ref, ...
                     M_linear, K, K2, cw, sensor_mask, reflection_order);
        p_time = p_time + reflection;  
    else         % reflection without correction
        reflection = Reflection1D(mgrid, medium, p_ref, ...
                     M_linear, K, cw, sensor_mask, reflection_order);
        p_time = p_time + reflection;  
    end
            
end


end