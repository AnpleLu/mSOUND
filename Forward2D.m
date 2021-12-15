function [p_time] = Forward2D(mgrid, medium, excit_ps, sensor_mask,...
                    reflection_order, varargin)
% DESCRIPTION:
% Simpson-like integration in the forward direction 

% USAGE:
% [p_time] = Forward2D(mgrid, medium, excit_ps, sensor_mask,...
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
% 'NRL'              apply non-reflecting layer at the boundary to reduce
%                    the spatial aliasing error
% 'animation'        visualize the simulation
% 'record_animation' generate a movie after the simulation

% OUTPUTS:
% p_time             waveforms recorded at the sensor positions

%%
% reshape the excitation signal for forward projection
excit_p = excitation(excit_ps, mgrid);
clear excit_ps

% define p_time to store the time-domain results
p_time = zeros(mgrid.num_t, sum(sum(sensor_mask)));

% if the first column of the sensor_mask is activated
if sum(sum(sensor_mask(:,1)))~=0
    % calculate the ending position of recorded signal
    sensor_mask_I2 = sum(sum(sensor_mask(:,1)));
    % save the time-domain signals
    p_time(:,1:sensor_mask_I2) = excit_p(:,sensor_mask(:,1)~=0);
end  

if reflection_order~=0
    % define p_ref to store the reflected part
    p_ref = zeros(mgrid.num_t,mgrid.num_x, mgrid.num_y+1);
end

% 1.05 is an empirical number, it can be higher if necessay
evanescent_fiter = 1.05;

% find medium properties
medium_cond = medium_case(medium);
MDMC = 0;
MDMA = 0;
MDM_ani = 0;
MDM_mov = 0;
if ~isempty(varargin)
    for input_index = 1:length(varargin)
        switch varargin{input_index}
            case 'correction'
            MDMC = 1;     % flag for correction
            case 'NRL'
            MDMA = 1;     % flag for non-reflecting layer
            
            case 'animation'
            MDM_ani = 1;  % flag for animation  
            p_ani = zeros(mgrid.num_t, mgrid.num_x, mgrid.num_y+1);
            p_ani(:,:,1) = excit_p;          
            figure
            imagesc (mgrid.y*1e3, mgrid.x*1e3, squeeze(p_ani(1,:,:)))
            xlabel ('x (mm)')
            ylabel ('y (mm)')
            axis image
            
            case 'record_animation'
            MDM_mov = 1;  % flag for recording the animation  
            MDM_ani = 1;
            p_ani = zeros(mgrid.num_t, mgrid.num_x, mgrid.num_y+1);   
            p_ani(:,:,1) = excit_p;  
            figure
            imagesc (mgrid.y*1e3, mgrid.x*1e3, squeeze(p_ani(1,:,:)))
            xlabel ('x (mm)')
            ylabel ('y (mm)')
            axis image    
            writerObj = VideoWriter('animation.avi');
            open(writerObj);            
        end
    end
end

%% absorption layer
if MDMA ==1
    M_gamma = repmat(medium.NRL_gamma./...
              (cosh(medium.NRL_alpha.*mgrid.abx_vec)).^2.*1i,mgrid.num_t,1).*...
              (mgrid.w'*ones(1,mgrid.num_x));
else 
    M_gamma = 0;
end  

%%
if medium_cond == 1
    
    disp( 'Linear homogeneous medium without dispersion')
    
    % calculate the M term on RHS of solution equation
    [~, M_linear, cw] = Mterm2D(mgrid, medium, medium_cond);
    
    % add non-reflecting layer
    M_linear = repmat(M_linear, 1, mgrid.num_x) + M_gamma; 
     
    % wave vector along the propagation direction
    Ktemp1 = mgrid.w'.^2*ones(1,mgrid.num_x)./medium.c0.^2 -...
             ones(mgrid.num_t,1)*mgrid.kx.^2 - ...
             M_linear; 
    K = sqrt(Ktemp1);
    clear Ktemp1
     
    if mod(mgrid.num_t,2) == 1
        K(mgrid.num_t/2+0.5:end,:) = -K(mgrid.num_t/2+0.5:end,:);
    else
        K(mgrid.num_t/2:end,:) = -K(mgrid.num_t/2:end,:);
    end
    K(K==0) = eps;        
        
    f = excit_p;
    clear  excit_p
     
    % Fourier transform of the excitation with regard of x and time t
    excit_F = fftshift(fft2(f)); 
     
    h = waitbar(0,'Forward projection, please wait...');
        
    for I = 2:mgrid.num_y+1
        
        Ktemp1 = (mgrid.w'./cw).^2*ones(1,mgrid.num_x) -...
                 evanescent_fiter*(ones(mgrid.num_t,1)*mgrid.kx.^2 + M_linear);    

        excit_F = excit_F.*exp(1i*K*mgrid.dy); % forward projection
        excit_F(isnan(excit_F)) = 0;  
        excit_F(real(Ktemp1)<=0) = 0;          % filter out evanescent wave
        f = real(ifft2(ifftshift(excit_F)));   % p0(z+mgrid.dy) 
         
        
        if max(max(max(abs(f))))==0
            disp ('The computation has lead to unphysical results. Please try reducing your step size.')
            break
        end          
        
        if MDM_ani % show the animation
            p_ani(:,:,I) = f;
            it = round(I*mgrid.dy/medium.c0/mgrid.dt);
            
            % in case a small computational time domain is applied
            it = max(mod(it, mgrid.num_t),1);
            
            imagesc (mgrid.y*1e3, mgrid.x*1e3, squeeze(p_ani(it,:,:)))
            axis image
            drawnow;  
            if MDM_mov % record the animation
                writeVideo(writerObj, getframe(gcf));
            end
        end
   
        % save the result according tod the sensor_mask
        if sum(sum(sensor_mask(:,I)))~=0
            
            % calculate the starting position of recorded signal
            sensor_mask_I1 = sum(sum(sensor_mask(:,1:I-1)))+1;
       
            % calculate the ending position of recorded signal
            sensor_mask_I2 = sensor_mask_I1-1 + ...
            sum(sum(sensor_mask(:,I)));
            % save the time-domain signals
            p_time(:,sensor_mask_I1:sensor_mask_I2) = f(:,sensor_mask(:,I)~=0);
        end  
         
        waitbar(I/mgrid.num_y)
    end
    close(h)

%%
elseif medium_cond == 2
    
    disp ('Nonlinear homogeneous medium')
      
    % calculate the M term on RHS of solution equation
    [M_nonlinear, M_linear, cw] = Mterm2D(mgrid, medium, medium_cond);      

    % add non-reflecting layer
    M_linear = repmat(M_linear, 1, mgrid.num_x) + M_gamma;
    
    % wave vector along the propagation direction
    Ktemp1 = mgrid.w'.^2*ones(1,mgrid.num_x)./medium.c0.^2 -...
             ones(mgrid.num_t,1)*mgrid.kx.^2 - ...
             M_linear;      
           
    K = sqrt(Ktemp1);
    clear Ktemp1

    if mod(mgrid.num_t,2) == 1
        K(mgrid.num_t/2+0.5:end,:) = -K(mgrid.num_t/2+0.5:end,:);
    else
        K(mgrid.num_t/2:end,:) = -K(mgrid.num_t/2:end,:);
    end
    K(K==0) = eps;

    % exponential term for the forward propagation
    expn2 = exp(0.5*1i*K*mgrid.dy); 
        
    f = excit_p;
    clear excit_p
    
    % Fourier transform of the excitation with regard of x and time
    excit_F = fftshift(fft2(f)); 
   
    
    h = waitbar(0,'Forward projection, please wait...');
        
    for I = 2:mgrid.num_y+1
        
        Ktemp1 = (mgrid.w'./cw).^2*ones(1,mgrid.num_x) -...
                  evanescent_fiter*ones(mgrid.num_t,1)*mgrid.kx.^2 - ...
                  evanescent_fiter*M_linear;                
         
        % Simpson 
        % 1   
        M = repmat(M_nonlinear, 1,mgrid.num_x).*fftshift(fft2(f.*f));
        F1 = excit_F.*expn2 + 0.5*mgrid.dy*M.*expn2./(2i.*K);  % P01(z+mgrid.dy/2)
        F1(isnan(F1)) = 0;
        F1(real(Ktemp1)<=0) = 0; 
    
        % 2
        f = real(ifft2(ifftshift(F1)));  % p0(z+mgrid.dy/2) 
        M1 = repmat(M_nonlinear, 1,mgrid.num_x).*fftshift(fft2(f.*f));
        F2 = excit_F.*expn2 + 0.25*mgrid.dy*expn2./(2i.*K).*(M + M1./expn2);% P12(z+mgrid.dy/2)
        clear M1  F1
        F2(isnan(F2)) = 0;
        F2(real(Ktemp1)<=0) = 0;   
  
        % 3   
        f  = real(ifft2(ifftshift(F2))); % p1(z+mgrid.dy/2)  
        M2 = repmat(M_nonlinear, 1,mgrid.num_x).*fftshift(fft2(f.*f));
        F3 = F2.*expn2 + 0.5*mgrid.dy*M2.*expn2./(2i.*K);  % P03(z+mgrid.dy)
        F3(isnan(F3)) = 0;
        F3(real(Ktemp1)<=0) = 0;  
    
        % 4     
        f  = real(ifft2(ifftshift(F3))); % p0(z+mgrid.dy)   
        clear  F3
        M3 = repmat(M_nonlinear, 1,mgrid.num_x).*fftshift(fft2(f.*f));
        F4 = F2.*expn2 + 0.25*mgrid.dy.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z+mgrid.dy)
        clear M3  F2 
        F4(isnan(F4)) = 0;
        F4(real(Ktemp1)<=0) = 0;    
      
        % 5   
        f  = real(ifft2(ifftshift(F4))); % p0(z+mgrid.dy) 
        M4 = repmat(M_nonlinear, 1,mgrid.num_x).*fftshift(fft2(f.*f));
        F5 = excit_F.*exp(1i*K*mgrid.dy) +...
             mgrid.dy/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*mgrid.dy)).*...
             exp(1i*K*mgrid.dy)./(2i.*K); % P25(z+mgrid.dy)
     
         clear M  M2  M4  F4
         F5(isnan(F5)) = 0;
         F5(real(Ktemp1)<=0) = 0; 
         excit_F = F5;
         f = real(ifft2(ifftshift(F5))); % p0(z+mgrid.dy) 
         clear F5
      
        if max(max(max(abs(f))))==0
            disp ('The computation has lead to unphysical results. Please try reducing your step size.')
            break
        end           
         
         if MDM_ani % show the animation
             p_ani(:,:,I) = f;
             it = round(I*mgrid.dy/medium.c0/mgrid.dt);
             
             % in case a small computational time domain is applied
             it = max(mod(it, mgrid.num_t),1);
             
             imagesc (mgrid.y*1e3, mgrid.x*1e3, squeeze(p_ani(it,:,:)))
             axis image
             drawnow;  
             if MDM_mov % record the animation
                 writeVideo(writerObj, getframe(gcf));
             end
         end
 %% save data   
         if sum(sum(sensor_mask(:,I)))~=0
             % calculate the starting position of recorded signal
             sensor_mask_I1 = sum(sum(sensor_mask(:,1:I-1)))+1;
       
             % calculate the ending position of recorded signal
             sensor_mask_I2 = sensor_mask_I1 -1 + ...
             sum(sum(sensor_mask(:,I)));
    
              % save the time-domain signals
              p_time(:,sensor_mask_I1:sensor_mask_I2) = f(:,sensor_mask(:,I)~=0);
         end
         waitbar(I/mgrid.num_y)
    end
    close(h)      
    
%%
elseif medium_cond == 3
    
    disp ('Inhomogeneous only in nonlinearity coefficient"')
        
    % normalized wave field
    f = excit_p;
    clear excit_p    
    
    % Fourier transform of the excitation with regard of x and time
    excit_F = fftshift(fft2(f)); 
        
    % calculate the M term on RHS of solution equation 
    % in this case, M_linear means lossy medium
    [M_nonlinear, M_linear, cw] = Mterm2D(mgrid, medium, medium_cond);
    
    % add anon-reflecting layer
    M_linear = repmat(M_linear, 1, mgrid.num_x) + M_gamma; 
    
    % calculate wavenumber along the propagation direction
    Ktemp1 = mgrid.w'.^2*ones(1,mgrid.num_x)./medium.c0.^2 -...
             ones(mgrid.num_t,1)*mgrid.kx.^2 -...
             M_linear; 
    K = sqrt(Ktemp1);
    clear Ktemp1
    if mod(mgrid.num_t,2) == 1
        K(mgrid.num_t/2+0.5:end,:) = -K(mgrid.num_t/2+0.5:end,:);
    else
        K(mgrid.num_t/2:end,:) = -K(mgrid.num_t/2:end,:);
    end
    K(K==0) = eps;            
        
    % exponential term for the forward propagation
    expn2 = exp(0.5*1i*K*mgrid.dy); 
    
    h = waitbar(0,'Forward projection, please wait...');

    for I = 2:mgrid.num_y+1
        
        Ktemp1 = mgrid.w'.^2*ones(1,mgrid.num_x)./medium.c0.^2 -...
                 evanescent_fiter*(ones(mgrid.num_t,1)*mgrid.kx.^2 + M_linear);    
             
        % Simpson 
        % 1   
        F2t = fftshift(fft(f.*f,[],1),1);
        % M(f(z))
        M = fftshift(fft(M_nonlinear(:,:,I-1).*F2t,[],2),2);
        F1 = excit_F.*expn2 + 0.5*mgrid.dy*M.*expn2./(2i.*K);  % P01(z+mgrid.dy/2)
        F1(isnan(F1)) = 0;
        F1(real(Ktemp1)<=0) = 0; 
    
        % 2
        f = real(ifft2(ifftshift(F1)));  % p0(z+mgrid.dy/2) 
        F2t = fftshift(fft(f.*f,[],1),1);  % fft(f^2)   
        % M(f0(z+mgrid.dy/2))  M1
        M1 = fftshift(fft(M_nonlinear(:,:,I-1).*F2t,[],2),2); 
        F2 = excit_F.*expn2 + 0.25*mgrid.dy*expn2./(2i.*K).*(M + M1./expn2);% P12(z+mgrid.dy/2)
        clear M1  F1
        F2(isnan(F2)) = 0;
        F2(real(Ktemp1)<=0) = 0;   
  
        % 3   
        f   = real(ifft2(ifftshift(F2))); % p1(z+mgrid.dy/2)  
        F2t = fftshift(fft(f.*f,[],1),1); % fft(f1^2)   
        % M(f1(z+mgrid.dy/2)) M2
        M2 = fftshift(fft(M_nonlinear(:,:,I-1).*F2t,[],2),2);
        F3 = F2.*expn2 + 0.5*mgrid.dy*M2.*expn2./(2i.*K);  % P03(z+mgrid.dy)
        F3(isnan(F3)) = 0;
        F3(real(Ktemp1)<=0) = 0;  
    
        % 4     
        f   = real(ifft2(ifftshift(F3))); % p0(z+mgrid.dy)   
        clear  F3
        F2t = fftshift(fft(f.*f,[],1),1); % fft(f2^2)
        % M(f0(z+mgrid.dy))  M3
        M3 = fftshift(fft(M_nonlinear(:,:,I).*F2t,[],2),2); 
        F4 = F2.*expn2 + 0.25*mgrid.dy.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z+mgrid.dy)
        clear M3  F2 
        F4(isnan(F4)) = 0;
        F4(real(Ktemp1)<=0) = 0;    
      
        % 5   
        f   = real(ifft2(ifftshift(F4))); % p0(z+mgrid.dy) 
        F2t = fftshift(fft(f.*f,[],1),1); % fft(f2^2)
        % M(f1(z+mgrid.dy)) M4
        M4 = fftshift(fft(M_nonlinear(:,:,I).*F2t,[],2),2);
        F5 = excit_F.*exp(1i*K*mgrid.dy) +...
             mgrid.dy/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*mgrid.dy)).*...
             exp(1i*K*mgrid.dy)./(2i.*K); % P25(z+mgrid.dy)
     
        clear M  M2  M4  F4
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
        excit_F = F5;
    
        f = real(ifft2(ifftshift(F5))); % p0(z+mgrid.dy) 
        clear F5
       
        if max(max(max(abs(f))))==0
            disp ('The computation has lead to unphysical results. Please try reducing your step size.')
            break
        end          
        
        if MDM_ani % show the animation
            p_ani(:,:,I) = f;
            it = round(I*mgrid.dy/medium.c0/mgrid.dt);
            
            % in case a small computational time domain is applied
            it = max(mod(it, mgrid.num_t),1);
            
            imagesc (mgrid.y*1e3, mgrid.x*1e3, squeeze(p_ani(it,:,:)))
            axis image
            drawnow;  
            if MDM_mov % record the animation
                writeVideo(writerObj, getframe(gcf));
            end
        end

        if sum(sum(sensor_mask(:,I)))~=0
            % calculate the starting position of recorded signal
            sensor_mask_I1 = sum(sum(sensor_mask(:,1:I-1)))+1;
       
            % calculate the ending position of recorded signal
            sensor_mask_I2 = sensor_mask_I1 -1 + ...
            sum(sum(sensor_mask(:,I)));
    
            % save the time-domain signals
            p_time(:,sensor_mask_I1:sensor_mask_I2) =...
                 f(:,sensor_mask(:,I)~=0);
        end  
       
        waitbar(I/mgrid.num_y)
    end
    close(h)
            
    
elseif medium_cond == 4 && MDMC == 0
    
    disp ('Linear medium without transmission corrections')
    
    % calculate the M term on RHS of solution equation
    [~, M_linear, cw] = Mterm2D(mgrid, medium, medium_cond);
    
    % add absoprtion layer
    M_linear = M_linear + repmat(M_gamma, 1,1,mgrid.num_y+1);
    
    % wavevector along the propagation direction
    Ktemp1 = (mgrid.w.'./medium.c0).^2*ones(1, mgrid.num_x) - ones(mgrid.num_t,1)*mgrid.kx.^2;  
    K  = sqrt(Ktemp1);
    clear Ktemp1
    if mod(mgrid.num_t,2) == 1
        K(mgrid.num_t/2+0.5:end,:) = -K(mgrid.num_t/2+0.5:end,:);
    else
        K(mgrid.num_t/2:end,:) = -K(mgrid.num_t/2:end,:);
    end    
    
    if reflection_order~=0
        Ktemp2 = zeros(mgrid.num_t, mgrid.num_x, mgrid.num_y+1);
        for it = 1:mgrid.num_t
            wt = mgrid.w(it)*ones(mgrid.num_x, mgrid.num_y+1);
            Ktemp2(it,:,:) = (wt./medium.c).^2 - mgrid.kx.'.^2*ones(1,mgrid.num_y+1);    
        end
        K2 = sqrt(Ktemp2);    
        clear  Ktemp2
    
        if mod(mgrid.num_t,2) == 1
            K2(mgrid.num_t/2+0.5:end,:) = -K2(mgrid.num_t/2+0.5:end,:);
        else
            K2(mgrid.num_t/2:end,:,:) = -K2(mgrid.num_t/2:end,:,:);
        end
        K2(K2==0) = eps;     
    else
        K2 = K;
    end    
    
    % exponential term for the forward propagation
    expn2 = exp(0.5*1i*K*mgrid.dy); 

    rho_rho = medium.rho;
    if length(medium.rho) == 1
        rho_rho = medium.rho*ones(mgrid.num_x, mgrid.num_y+1);
    end        
    
    % normalized wave field
    f = excit_p./(ones(mgrid.num_t,1)*sqrt(rho_rho(:,1).'));
        
    % Fourier transform of the excitation with regard of x and time
    excit_F = fftshift(fft2(f)); 

    h = waitbar(0,'Forward projection, please wait...');
    tic    
    for I = 2:mgrid.num_y+1
           
        Ktemp1 = (mgrid.w'*ones(1,mgrid.num_x)./cw(:,:,I)).^2 - ...
                  evanescent_fiter*ones(mgrid.num_t,1)*mgrid.kx.^2;    
         
        %Simpson 
        % 1   
        Ft  = fftshift(fft(f,[],1),1);
        % M(f(z))
        M = fftshift(fft((M_linear(:,:,I-1)).*Ft,[],2),2);
        F1 = excit_F.*expn2 + 0.5*mgrid.dy*M.*expn2./(2i.*K);  % P01(z+mgrid.dy/2)
        F1(isnan(F1)) = 0;
        F1(real(Ktemp1)<=0) = 0; 
    
        % 2
        f = real(ifft2(ifftshift(F1)));  % p0(z+mgrid.dy/2) 
        Ft  = fftshift(fft(f,[],1),1);     % fft(f) 
        % M(f0(z+mgrid.dy/2))  M1
        M1 = fftshift(fft((M_linear(:,:,I-1)).*Ft,[],2),2);
        F2 = excit_F.*expn2 + 0.25*mgrid.dy*expn2./(2i.*K).*(M + M1./expn2);% P12(z+mgrid.dy/2)
        clear M1  F1
        F2(isnan(F2)) = 0;
        F2(real(Ktemp1)<=0) = 0;   
  
        % 3   
        f   = real(ifft2(ifftshift(F2))); % p1(z+mgrid.dy/2)  
        Ft  = fftshift(fft(f,[],1),1);    % fft(f1) 
        % M(f1(z+mgrid.dy/2)) M2
        M2 = fftshift(fft((M_linear(:,:,I-1)).*Ft,[],2),2);
        F3 = F2.*expn2 + 0.5*mgrid.dy*M2.*expn2./(2i.*K);  % P03(z+mgrid.dy)
        F3(isnan(F3)) = 0;
        F3(real(Ktemp1)<=0) = 0;  
    
        % 4     
        f   = real(ifft2(ifftshift(F3))); % p0(z+mgrid.dy)   
        clear  F3
        Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 

        % M(f0(z+mgrid.dy))  M3
        M3 = fftshift(fft((M_linear(:,:,I)).*Ft,[],2),2);
        F4 = F2.*expn2 + 0.25*mgrid.dy.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z+mgrid.dy)
        clear M3  F2 
        F4(isnan(F4)) = 0;
        F4(real(Ktemp1)<=0) = 0;    
      
        % 5   
        f   = real(ifft2(ifftshift(F4))); % p0(z+mgrid.dy) 
        Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
        % M(f1(z+mgrid.dy)) M4
        M4 = fftshift(fft((M_linear(:,:,I)).*Ft,[],2),2);
        F5 = excit_F.*exp(1i*K*mgrid.dy) +...
             mgrid.dy/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*mgrid.dy)).*...
             exp(1i*K*mgrid.dy)./(2i.*K); % P25(z+mgrid.dy)
     
        clear M  M2  M4  F4
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
%%
        rho2 = ones(mgrid.num_t,1)*rho_rho(:,I).'; 
        
        if reflection_order~=0
            rho1 = (repmat(rho_rho(:,I-1)', mgrid.num_t,1));
            T_rhoc = 2.*rho2.*cw(:,:,I)./(rho1.*cw(:,:,I-1) + rho2.*cw(:,:,I));  % transmission coefficient
            P_in = excit_F.*exp(1i.*K2(:,:,I)*mgrid.dy);
            p_ref(:,:,I-1) = real(ifft2(ifftshift(P_in))).*(T_rhoc-1);  
            clear rho1  T_rhoc P_in
        end
            
        excit_F = F5;
        f = real(ifft2(ifftshift(F5))); % p0(z+mgrid.dy) 
        clear F5
    
        if max(max(max(abs(f))))==0
            disp ('The computation has lead to unphysical results. Please try reducing your step size.')
            break
        end          
        
        % recover the normalized wave field
        p_I = sqrt(rho2).*f;  
        clear rho2
  
        if MDM_ani % show the animation
            p_ani(:,:,I) = p_I;
            it = round(I*mgrid.dy/medium.c0/mgrid.dt);
            
            % in case a small computational time domain is applied
            it = max(mod(it, mgrid.num_t),1);
            
            imagesc (mgrid.y*1e3, mgrid.x*1e3, squeeze(p_ani(it,:,:)))
            axis image
            drawnow;  
            if MDM_mov % record the animation
                writeVideo(writerObj, getframe(gcf));
            end
        end         

        if sum(sum(sensor_mask(:,I)))~=0
            % calculate the starting position of recorded signal
            sensor_mask_I1 = sum(sum(sensor_mask(:,1:I-1)))+1;
       
            % calculate the ending position of recorded signal
            sensor_mask_I2 = sensor_mask_I1 -1 + ...
            sum(sum(sensor_mask(:,I)));
    
            % save the time-domain signals
            p_time(:,sensor_mask_I1:sensor_mask_I2) = p_I(:,sensor_mask(:,I)~=0);
        end        
        
        clear p_I
        waitbar(I/mgrid.num_y)
    end
    close(h)    
    toc
%%
elseif medium_cond ==4 && MDMC == 1 % add correction
    
    disp ('Linear medium with transmission corrections')
    
    % calculate the M term on RHS of solution equation
    [~, M_linear, Ktemp1, Ktemp2, cw] = Mterm2D_MMDM(mgrid, medium, medium_cond);
    
    % add absorption layer
    M_linear = M_linear + repmat(M_gamma, 1,1,mgrid.num_y+1);
    
    K  = sqrt(Ktemp1);
    K2 = sqrt(Ktemp2);
    
    if mod(mgrid.num_t,2) == 1
        K(mgrid.num_t/2+0.5:end,:)  = -K(mgrid.num_t/2+0.5:end,:);
        K2(mgrid.num_t/2+0.5:end,:) = -K2(mgrid.num_t/2+0.5:end,:);
    else
        K(mgrid.num_t/2:end,:)  = -K(mgrid.num_t/2:end,:);
        K2(mgrid.num_t/2:end,:,:) = -K2(mgrid.num_t/2:end,:,:);
    end
    K(K==0)   = eps;
    K2(K2==0) = eps;
    
    f = excit_p;
    clear excit_p
    
    % Fourier transform of the excitation with regard of x and time
    excit_F = fftshift(fft2(f)); 

    % exponential term for the forward propagation
    expn2 = exp(0.5*1i*K*mgrid.dy); 
    
    rho_rho = medium.rho;
    if length(medium.rho) == 1
        rho_rho = medium.rho*ones(mgrid.num_x, mgrid.num_y+1);
    end      
    
    h = waitbar(0,'Forward projection, please wait...');
    tic
    for I = 2:mgrid.num_y+1

        Ktemp = (mgrid.w'*ones(1,mgrid.num_x)./cw(:,:,I)).^2 -...
                evanescent_fiter*ones(mgrid.num_t,1)*mgrid.kx.^2;    
        % Simpson 

        % 1   
        Ft  = fftshift(fft(f,[],1),1);
        % M(f(z))
        M1 = fftshift(fft((M_linear(:,:,I-1)).*Ft,[],2),2);
        F1 = excit_F.*expn2 + 0.5*mgrid.dy*M1.*expn2./(2i.*K);  % P01(z+dz/2)
        F1(isnan(F1)) = 0;
        F1(real(Ktemp)<=0) = 0;   
    
        % 2
        f = real(ifft2(ifftshift(F1)));  % p0(z+dz/2) 
        clear  F1
        Ft  = fftshift(fft(f,[],1),1);     % fft(f) 
        % M(f0(z+dz/2))  M1
        M2 = fftshift(fft((M_linear(:,:,I-1)).*Ft,[],2),2);
        F2 = excit_F.*expn2 + 0.25*mgrid.dy*expn2./(2i.*K).*(M1 + M2./expn2);% P12(z+dz/2)
        clear M2
        F2(isnan(F2)) = 0;
        F2(real(Ktemp)<=0) = 0;  
  
        % 3   
        f   = real(ifft2(ifftshift(F2))); % p1(z+dz/2)  
        Ft  = fftshift(fft(f,[],1),1);    % fft(f1) 
        % M(f1(z+dz/2)) M2
        M3 = fftshift(fft((M_linear(:,:,I-1)).*Ft,[],2),2);
        F3 = F2.*expn2 + 0.5*mgrid.dy*M3.*expn2./(2i.*K);  % P03(z+dz)
        F3(isnan(F3)) = 0;
        F3(real(Ktemp)<=0) = 0;  
    
        % 4     
        f   = real(ifft2(ifftshift(F3))); % p0(z+dz)   
        clear F3 
        Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
        % M(f0(z+dz))  M3
        M4 = fftshift(fft((M_linear(:,:,I)).*Ft,[],2),2);
        F4 = F2.*expn2 + 0.25*mgrid.dy.*(M3 + M4./expn2).*expn2./(2i.*K); % P14(z+dz)
        F4(isnan(F4)) = 0;
        F4(real(Ktemp)<=0) = 0;    
        
        % 5   
        f   = real(ifft2(ifftshift(F4))); % p0(z+dz) 
        clear F4  
        Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
        % M(f1(z+dz)) M4
        M5 = fftshift(fft((M_linear(:,:,I)).*Ft,[],2),2);
        F5 = excit_F.*exp(1i*K*mgrid.dy) +...
             mgrid.dy/6.0.*(M1 + 4*M3./expn2 + M5./exp(1i*K*mgrid.dy)).*...
             exp(1i*K*mgrid.dy)./(2i.*K); % P25(z+mgrid.dy)
        clear M1  M3 M5
        F5(isnan(F5)) = 0;
        F5(real(Ktemp)<=0) = 0; 
        %% phase 
        c1 = cw(:,:,I-1);
        c2 = cw(:,:,I); 
        rho1 = ones(mgrid.num_t,1)*rho_rho(:,I-1).'; 
        rho2 = ones(mgrid.num_t,1)*rho_rho(:,I).'; 

        wx2 = (mgrid.w.*mgrid.w).'*ones(1,mgrid.num_x);
        e1  = exp(1i*(K - K2(:,:,I))*mgrid.dy);
        k_dif2 = (wx2./medium.c0./medium.c0 - wx2./c2./c2)*mgrid.dy;      
        k_dif  = (wx2./medium.c0./medium.c0 - wx2./c1./c1)*mgrid.dy; 
        clear  wx2  
        correction = excit_F.*exp(1i*K2(:,:,I)*mgrid.dy).*...
                     (1-e1.*(1+k_dif./4i./K)./(1-k_dif2./4i./K));
       
        F5 = F5 + correction;       
        clear correction k_dif2  k_dif  el

        F5(isnan(F5)) = 0;
        F5(real(Ktemp)<=0) = 0; 
        %% amplitude correction    
        f = real(ifft2(ifftshift(F5)));
        T_rhoc = 2.*rho2.*c2./(rho1.*c1 + rho2.*c2);  % transmission coefficient
        f = f.*T_rhoc; 

        F5 = fftshift(fft2(f));
%%    
        if reflection_order~=0
            P_in = excit_F.*exp(1i.*K2(:,:,I)*mgrid.dy);
            p_ref(:,:,I-1) = real(ifft2(ifftshift(P_in))).*(T_rhoc-1);    
            clear P_in
        end
        %%
        excit_F = F5;
        f = real(ifft2(ifftshift(F5))); % p0(z+mgrid.dy) 

        if max(max(max(abs(f))))==0
            disp ('The computation has lead to unphysical results. Please try reducing your step size.')
            break
        end          
        
        
        clear F5
        
        if MDM_ani % show the animation
            p_ani(:,:,I) = f;
            it = round(I*mgrid.dy/medium.c0/mgrid.dt);
            
            % in case a small computational time domain is applied
            it = max(mod(it, mgrid.num_t),1);
            
            imagesc (mgrid.y*1e3, mgrid.x*1e3, squeeze(p_ani(it,:,:)))
            axis image
            drawnow;  
            if MDM_mov % record the animation
                writeVideo(writerObj, getframe(gcf));
            end
        end           

        if sum(sum(sensor_mask(:,I)))~=0
           % calculate the starting position of recorded signal
           sensor_mask_I1 = sum(sum(sensor_mask(:,1:I-1)))+1;
       
           % calculate the ending position of recorded signal
           sensor_mask_I2 = sensor_mask_I1 -1 + ...
           sum(sum(sensor_mask(:,I)));
    
           % save the time-domain signals
           p_time(:,sensor_mask_I1:sensor_mask_I2) = f(:,sensor_mask(:,I)~=0);
        end
       
        waitbar(I/mgrid.num_y)
    
    end
    toc
    close(h)    
    
   
 %%        
elseif medium_cond == 5 && MDMC == 0
    
%     disp ('Nonlinear media without transimission corrections')
    % calculate the M term on RHS of solution equation
    [M_nonlinear, M_linear,cw] = Mterm2D(mgrid, medium, medium_cond); 
    
    % add non-reflecting layer
    M_linear = M_linear + repmat(M_gamma, 1,1,mgrid.num_y+1);
      
    % calculate the wavevector in propagation direction
    Ktemp1 =(mgrid.w.'./medium.c0).^2*ones(1, mgrid.num_x) - ones(mgrid.num_t,1)*mgrid.kx.^2;  
    K  = sqrt(Ktemp1);
    clear Ktemp1
    
    if mod(mgrid.num_t,2) == 1
        K(mgrid.num_t/2+0.5:end,:)  = -K(mgrid.num_t/2+0.5:end,:);
    else
        K(mgrid.num_t/2:end,:)  = -K(mgrid.num_t/2:end,:);
    end
    K(K==0) = eps;     

    
    if reflection_order~=0
        Ktemp2 = zeros(mgrid.num_t, mgrid.num_x, mgrid.num_y+1);
        for it = 1:mgrid.num_t
            wt = mgrid.w(it)*ones(mgrid.num_x, mgrid.num_y+1);
            Ktemp2(it,:,:) = (wt./medium.c).^2 - mgrid.kx.'.^2*ones(1,mgrid.num_y+1);    
        end
        K2 = sqrt(Ktemp2);    
        clear  Ktemp2
    
        if mod(mgrid.num_t,2) == 1
            K2(mgrid.num_t/2+0.5:end,:) = -K2(mgrid.num_t/2+0.5:end,:);
        else
            K2(mgrid.num_t/2:end,:,:) = -K2(mgrid.num_t/2:end,:,:);
        end
        K2(K2==0) = eps;     
    else
        K2 = K;
    end      

    % exponential term for the forward propagation
    expn2 = exp(0.5*1i*K*mgrid.dy); 
 
    rho_rho = medium.rho;
    if length(medium.rho) == 1
        rho_rho = medium.rho*ones(mgrid.num_x, mgrid.num_y+1);
    end    
    
    % normalized wave field
    f = excit_p./(ones(mgrid.num_t,1)*sqrt(rho_rho(:,1).'));
        
    % Fourier transform of the excitation with regard of x and time
    excit_F = fftshift(fft2(f)); 
        
    h = waitbar(0,'Forward projection, please wait...');
        
    for I = 2:mgrid.num_y+1
        
        Ktemp1 = mgrid.w'.^2*ones(1,mgrid.num_x)./cw(:,:,I).^2 -...
                 evanescent_fiter*ones(mgrid.num_t,1)*mgrid.kx.^2;    
         
        % Simpson 
        % 1   
        Ft  = fftshift(fft(f,[],1),1);
        F2t = fftshift(fft(f.*f,[],1),1);
        % M(f(z))
        M = fftshift(fft((M_linear(:,:,I-1)).*Ft,[],2),2);
        M = M + fftshift(fft(M_nonlinear(:,:,I-1).*F2t,[],2),2);
        F1 = excit_F.*expn2 + 0.5*mgrid.dy*M.*expn2./(2i.*K);  % P01(z+mgrid.dy/2)
        F1(isnan(F1)) = 0;
        F1(real(Ktemp1)<=0) = 0; 
    
        % 2
        f = real(ifft2(ifftshift(F1)));  % p0(z+mgrid.dy/2) 
        Ft  = fftshift(fft(f,[],1),1);     % fft(f) 
        F2t = fftshift(fft(f.*f,[],1),1);  % fft(f^2)   
        % M(f0(z+mgrid.dy/2))  M1
        M1 = fftshift(fft((M_linear(:,:,I-1)).*Ft,[],2),2);
        M1 = M1 +  fftshift(fft(M_nonlinear(:,:,I-1).*F2t,[],2),2); 
        F2 = excit_F.*expn2 + 0.25*mgrid.dy*expn2./(2i.*K).*(M + M1./expn2);% P12(z+mgrid.dy/2)
        clear M1  F1
        F2(isnan(F2)) = 0;
        F2(real(Ktemp1)<=0) = 0;   
  
        % 3   
        f   = real(ifft2(ifftshift(F2))); % p1(z+mgrid.dy/2)  
        Ft  = fftshift(fft(f,[],1),1);    % fft(f1) 
        F2t = fftshift(fft(f.*f,[],1),1); % fft(f1^2)   
        % M(f1(z+mgrid.dy/2)) M2
        M2 = fftshift(fft((M_linear(:,:,I-1)).*Ft,[],2),2);
        M2 = M2 + fftshift(fft(M_nonlinear(:,:,I-1).*F2t,[],2),2);
        F3 = F2.*expn2 + 0.5*mgrid.dy*M2.*expn2./(2i.*K);  % P03(z+mgrid.dy)
        F3(isnan(F3)) = 0;
        F3(real(Ktemp1)<=0) = 0;  
    
        % 4     
        f   = real(ifft2(ifftshift(F3))); % p0(z+mgrid.dy)   
        clear  F3
        Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
        F2t = fftshift(fft(f.*f,[],1),1); % fft(f2^2)
        % M(f0(z+mgrid.dy))  M3
        M3 = fftshift(fft((M_linear(:,:,I)).*Ft,[],2),2);
        M3 = M3 + fftshift(fft(M_nonlinear(:,:,I).*F2t,[],2),2); 
        F4 = F2.*expn2 + 0.25*mgrid.dy.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z+mgrid.dy)
        clear M3  F2 
        F4(isnan(F4)) = 0;
        F4(real(Ktemp1)<=0) = 0;    
      
        % 5   
        f   = real(ifft2(ifftshift(F4))); % p0(z+mgrid.dy) 
        Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
        F2t = fftshift(fft(f.*f,[],1),1); % fft(f2^2)
        % M(f1(z+mgrid.dy)) M4
        M4 = fftshift(fft((M_linear(:,:,I)).*Ft,[],2),2);
        M4 = M4 + fftshift(fft(M_nonlinear(:,:,I).*F2t,[],2),2);
        F5 = excit_F.*exp(1i*K*mgrid.dy) +...
             mgrid.dy/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*mgrid.dy)).*...
             exp(1i*K*mgrid.dy)./(2i.*K); % P25(z+mgrid.dy)
     
        clear M  M2  M4  F4
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
%%  
        if reflection_order~=0
            c1 = cw(:,:,I-1);
            c2 = cw(:,:,I); 
            rho1 = ones(mgrid.num_t,1)*rho_rho(:,I-1).'; 
            rho2 = ones(mgrid.num_t,1)*rho_rho(:,I).';             

            T_rhoc = 2.*rho2.*c2./(rho1.*c1 + rho2.*c2);   % transmission coefficient
            P_in = excit_F.*exp(1i.*K2(:,:,I)*mgrid.dy);              
            p_ref(:,:,I-1) = real(ifft2(ifftshift(P_in))).*(T_rhoc-1);
        end
                        
        excit_F = F5;
        f = real(ifft2(ifftshift(F5))); % p0(z+mgrid.dy) 
        clear F5
        
        if max(max(max(abs(f))))==0
            disp ('The computation has lead to unphysical results. Please try reducing your step size.')
            break
        end          
        
        % recover the normalized wave field
        p_I = (ones(mgrid.num_t,1)*sqrt(rho_rho(:,I).')).*f;    
         
        if MDM_ani % show the animation
            p_ani(:,:,I) = p_I;
            
            it = round(I*mgrid.dy/min(cw(1,:,I))/mgrid.dt);
            
            % in case a small computational time domain is applied
            it = max(mod(it, mgrid.num_t),1);
            
            imagesc (mgrid.y*1e3, mgrid.x*1e3, squeeze(p_ani(it,:,:)))
            axis image
            drawnow;  
            if MDM_mov % record the animation
                writeVideo(writerObj, getframe(gcf));
            end
        end          
        
        if sum(sum(sensor_mask(:,I)))~=0
            % calculate the starting position of recorded signal
            sensor_mask_I1 = sum(sum(sensor_mask(:,1:I-1)))+1;
       
            % calculate the ending position of recorded signal
            sensor_mask_I2 = sensor_mask_I1 -1 + ...
            sum(sum(sensor_mask(:,I)));
    
            % save the time-domain signals
            p_time(:,sensor_mask_I1:sensor_mask_I2) = p_I(:,sensor_mask(:,I)~=0);
        
        end
        clear p_I
        waitbar(I/mgrid.num_y)
     end
     close(h)

%%
elseif medium_cond == 5 && MDMC == 1
    
%     disp ('Nonlinear media with transmission corrections')
    
    % calculate the M term on RHS of solution equation
    [M_nonlinear, M_linear, Ktemp1, Ktemp2, cw] = Mterm2D_MMDM(mgrid, medium, medium_cond);
    
    % add absoprtion layer
    M_linear = M_linear + repmat(M_gamma, 1,1,mgrid.num_y+1);
    
    % calculate wave number along the propagation direction
    K  = sqrt(Ktemp1);
    clear Ktemp1
    K2 = sqrt(Ktemp2);
    clear Ktemp2
    
    if mod(mgrid.num_t,2) == 1
        K(mgrid.num_t/2+0.5:end,:)  = -K(mgrid.num_t/2+0.5:end,:);
        K2(mgrid.num_t/2+0.5:end,:) = -K2(mgrid.num_t/2+0.5:end,:);
    else
        K(mgrid.num_t/2:end,:)  = -K(mgrid.num_t/2:end,:);
        K2(mgrid.num_t/2:end,:,:) = -K2(mgrid.num_t/2:end,:,:);
    end
    K(K==0)   = eps;
    K2(K2==0) = eps;
    
    f = excit_p;
    clear excit_p
    
    % Fourier transform of the excitation with regard of x and time
    excit_F = fftshift(fft2(f)); 

    % exponential term for the forward propagation
    expn2 = exp(0.5*1i*K*mgrid.dy); 
    
    rho_rho = medium.rho;
    if length(medium.rho) == 1
        rho_rho = medium.rho*ones(mgrid.num_x, mgrid.num_y+1);
    end     
    
    h = waitbar(0,'Forward projection, please wait...');
    
    for I = 2:mgrid.num_y+1
        
        Ktemp = mgrid.w'.^2*ones(1,mgrid.num_x)./cw(:,:,I).^2 -...
                evanescent_fiter*ones(mgrid.num_t,1)*mgrid.kx.^2;    
        % Simpson 

        % 1   
        Ft  = fftshift(fft(f,[],1),1);
        F2t = fftshift(fft(f.*f,[],1),1);
        % M(f(z))
        M1 = fftshift(fft((M_linear(:,:,I-1)).*Ft,[],2),2);
        M1 = M1 + fftshift(fft(M_nonlinear(:,:,I-1).*F2t,[],2),2);
        F1 = excit_F.*expn2 + 0.5*mgrid.dy*M1.*expn2./(2i.*K);  % P01(z+dz/2)
        F1(isnan(F1)) = 0;
        F1(real(Ktemp)<=0) = 0; 
    
        % 2
        f = real(ifft2(ifftshift(F1)));  % p0(z+dz/2) 
        clear  F1
        Ft  = fftshift(fft(f,[],1),1);     % fft(f) 
        F2t = fftshift(fft(f.*f,[],1),1);  % fft(f^2)   
        % M(f0(z+dz/2))  M1
        M2 = fftshift(fft((M_linear(:,:,I-1)).*Ft,[],2),2);
        M2 = M2 + fftshift(fft(M_nonlinear(:,:,I-1).*F2t,[],2),2); 
        F2 = excit_F.*expn2 + 0.25*mgrid.dy*expn2./(2i.*K).*(M1 + M2./expn2);% P12(z+dz/2)
        clear M2
        F2(isnan(F2)) = 0;
        F2(real(Ktemp)<=0) = 0;  
  
        % 3   
        f   = real(ifft2(ifftshift(F2))); % p1(z+dz/2)  
        Ft  = fftshift(fft(f,[],1),1);    % fft(f1) 
        F2t = fftshift(fft(f.*f,[],1),1); % fft(f1^2)   
        % M(f1(z+dz/2)) M2
        M3 = fftshift(fft((M_linear(:,:,I-1)).*Ft,[],2),2);
        M3 = M3 + fftshift(fft(M_nonlinear(:,:,I-1).*F2t,[],2),2);
        F3 = F2.*expn2 + 0.5*mgrid.dy*M3.*expn2./(2i.*K);  % P03(z+dz)
        F3(isnan(F3)) = 0;
        F3(real(Ktemp)<=0) = 0;  
    
        % 4     
        f   = real(ifft2(ifftshift(F3))); % p0(z+dz)   
        clear F3 
        Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
        F2t = fftshift(fft(f.*f,[],1),1); % fft(f2^2)
        % M(f0(z+dz))  M3
        M4 = fftshift(fft((M_linear(:,:,I)).*Ft,[],2),2);
        M4 = M4 + fftshift(fft(M_nonlinear(:,:,I).*F2t,[],2),2); 
        F4 = F2.*expn2 + 0.25*mgrid.dy.*(M3 + M4./expn2).*expn2./(2i.*K); % P14(z+dz)
        F4(isnan(F4)) = 0;
        F4(real(Ktemp)<=0) = 0;    
      
        % 5   
        f   = real(ifft2(ifftshift(F4))); % p0(z+dz) 
        clear F4  F2
        Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
        F2t = fftshift(fft(f.*f,[],1),1); % fft(f2^2)
        % M(f1(z+dz)) M4
        M5 = fftshift(fft((M_linear(:,:,I)).*Ft,[],2),2);
        M5 = M5 + fftshift(fft(M_nonlinear(:,:,I).*F2t,[],2),2);
        F5 = excit_F.*exp(1i*K*mgrid.dy) +...
             mgrid.dy/6.0.*(M1 + 4*M3./expn2 + M5./exp(1i*K*mgrid.dy)).*...
             exp(1i*K*mgrid.dy)./(2i.*K); % P25(z+mgrid.dy)
        clear M1  M3 M5
        F5(isnan(F5)) = 0;
        F5(real(Ktemp)<=0) = 0;  
        %% phase 
        c1 = cw(:,:,I-1);
        c2 = cw(:,:,I);    
        rho1 = ones(mgrid.num_t,1)*rho_rho(:,I-1).'; 
        rho2 = ones(mgrid.num_t,1)*rho_rho(:,I).';      
        
        wx2 = (mgrid.w.*mgrid.w).'*ones(1,mgrid.num_x);
        e1 = exp(1i*(K - K2(:,:,I))*mgrid.dy);
        k_dif2 = (wx2./medium.c0./medium.c0 - wx2./c2./c2)*mgrid.dy;      
        k_dif  = (wx2./medium.c0./medium.c0 - wx2./c1./c1)*mgrid.dy; 
        clear  wx2  
        correction = excit_F.*exp(1i*K2(:,:,I)*mgrid.dy).*...
                     (1-e1.*(1+k_dif./4i./K)./(1-k_dif2./4i./K));

        F5 = F5 + correction;  
        clear correction k_dif2  k_dif  el

        F5(isnan(F5)) = 0;
        F5(real(Ktemp)<=0) = 0; 
        %% amplitude correction    
        f = real(ifft2(ifftshift(F5)));
        
        T_rhoc = 2.*rho2.*c2./(rho1.*c1 + rho2.*c2);  %%layered medium
        f = f.*T_rhoc;
        F5 = fftshift(fft2(f));
        
%%  add reflection  
        if reflection_order~=0
            P_in = excit_F.*exp(1i.*K2(:,:,I)*mgrid.dy);
            p_ref(:,:,I-1) = real(ifft2(ifftshift(P_in))).*(T_rhoc-1);   
        end
          
%%
        excit_F = F5;
        f = real(ifft2(ifftshift(F5))); % p0(z+mgrid.dy) 

        if max(max(max(abs(f))))==0
            disp ('The computation has lead to unphysical results. Please try reducing your step size.')
            break
        end          
        
        if MDM_ani % show the animation
            p_ani(:,:,I) = f;
            it = round(I*mgrid.dy/medium.c0/mgrid.dt);
            
            % in case a small computational time domain is applied
            it = max(mod(it, mgrid.num_t),1);
            
            imagesc (mgrid.y*1e3, mgrid.x*1e3, squeeze(p_ani(it,:,:)))
            axis image
            drawnow;  
            if MDM_mov % record the animation
                writeVideo(writerObj, getframe(gcf));
            end
        end          

       if sum(sum(sensor_mask(:,I)))~=0
           % calculate the starting position of recorded signal
           sensor_mask_I1 = sum(sum(sensor_mask(:,1:I-1)))+1;
       
           % calculate the ending position of recorded signal
           sensor_mask_I2 = sensor_mask_I1 -1 + ...
           sum(sum(sensor_mask(:,I)));
    
           % save the time-domain signals
           p_time(:,sensor_mask_I1:sensor_mask_I2) = f(:,sensor_mask(:,I)~=0);
       end  
        
        waitbar(I/mgrid.num_y)
    end
    close(h)
    
end

% save the movie frames to disk
if MDM_mov
    close(writerObj);
end

% flag for reflection
if (max(max(medium.c))   ~= min(min(medium.c))) || ...
   (max(max(medium.rho)) ~= min(min(medium.rho)))

    reflection_option = 1; 
else
    reflection_option = 0;
end
    
% calcuate the reflections
if reflection_order ~=0 && (reflection_option ==1)
    if MDMC == 1 % reflection with correction

        reflection = MReflection2D(mgrid, medium, p_ref, ...
                     M_linear, K, K2, cw, sensor_mask, reflection_order);
        p_time = p_time + reflection;  
    else  % reflection without correction
        reflection = Reflection2D(mgrid, medium, p_ref, ...
                     M_linear, K, cw, sensor_mask, reflection_order);
        p_time = p_time + reflection;  
    end
            
end


end