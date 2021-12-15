function [p_time] = Backward3D(mgrid, medium, excit_ps, sensor_mask,...
                    reflection_order, varargin)
% DESCRIPTION:
% Simpson-like intergration in the backward direction 

% USAGE:
% [p_time] = Backward3D(mgrid, medium, excit_ps, sensor_mask,...
%                       reflection_order, varargin)

% INPUTS:
% mgrid              structure to define the computational domain
% medium             medium properties
% excit_ps           excitation signals
% sensor_mask        a set of Cartesian points where the pressure is recorded

% OPTIONAL INOUTS
% 'correction'       apply phase and amplitude corrections for strongly
%                    heterogeneous media
% 'NRL'              apply non-reflecting layer at the boundary to reduce
% the spatial aliasing error
% 'no_refl_correction' turnoff the reflection correction in backward
% projection

% OUTPUTS:
% p_time             waveforms recorded at the sensor positions

%%
% reshape the excitation signal for forward projection
excit_p = excitation(excit_ps, mgrid);

% 1.05 is an empirical number, it can be higher if necessay
evanescent_fiter = 1.05;

% define p_total to store the time-domain results
p_time = zeros(mgrid.num_t, sum(sum(sum(sensor_mask))));

% if the first column of the sensor_mask is activated
if sum(sum(sensor_mask(:,:,end)))~=0
    % calculate the ending position of recorded signal
    sensor_mask_I2 = sum(sum(sensor_mask(:,:,end)));
    % save the time-domain signals
    p_time(:,end-sensor_mask_I2+1:end) = excit_p(:,sensor_mask(:,:,end)~=0);
end  

if reflection_order~=0
    % define p_ref to store the reflected part
    p_ref = zeros(mgrid.num_t, mgrid.num_x,mgrid.num_y, mgrid.num_z+1);
end

% find medium properties
medium_cond = medium_case(medium);
MDMC = 0; % flag for correction
MDMA = 0; % flag for non-rflecting layer
MDM_ani = 0;
MDM_mov = 0;
MDMS = 1;
if ~isempty(varargin)
    for input_index = 1:length(varargin)
        switch varargin{input_index}
            case 'correction'
            MDMC = 1;     % flag for correction
            case 'NRL'
            MDMA = 1;     % flag non-rflecting layer
            case 'animation'
            MDM_ani = 1;  % flag for animation  
            p_ani = zeros(mgrid.num_t, mgrid.num_x, mgrid.num_z+1);
            p_ani(:,:,end) = excit_p(:,:,round(mgrid.num_y/2));
            figure
            imagesc (mgrid.z*1e3, mgrid.x*1e3, squeeze(p_ani(end,:,:)))
            xlabel ('x (mm)')
            ylabel ('y (mm)')
            
            case 'record_animation'
            MDM_mov = 1;  % flag for recording the animation  
            MDM_ani = 1;
            figure
            xlabel ('x (mm)')
            ylabel ('y (mm)')
            p_ani = zeros(mgrid.num_t, mgrid.num_x, mgrid.num_z+1); 
            p_ani(:,:,end) = excit_p(:,:,round(mgrid.num_y/2));
            figure
            imagesc (mgrid.z*1e3, mgrid.x*1e3, squeeze(p_ani(end,:,:)))            
            writerObj = VideoWriter('animation.avi');
            open(writerObj); 
            
            case 'no_refl_correction'  
            MDMS = 0; % turn off the reflection correction             
        end
    end
end

% calculate absorption layer properties
if MDMA ==1
    M_gamma = zeros(1, mgrid.num_x, mgrid.num_y);
    M_gamma(1,:,:) = repmat(medium.NRL_gamma./...
                     (cosh(medium.NRL_alpha.*mgrid.abx_vec.')).^2,1,mgrid.num_y) + ...
                     repmat(medium.NRL_gamma./...
                     (cosh(medium.NRL_alpha.*mgrid.aby_vec)).^2,mgrid.num_x,1);
    M_gamma = -1i*repmat(M_gamma,mgrid.num_t, 1,1).*...
              repmat(mgrid.w',1, mgrid.num_x, mgrid.num_y);                    
else  
    M_gamma = 0;
end  

backward_dz = -mgrid.dz;

if medium_cond == 1
    
    disp( 'Linear homogeneous medium')
     
    % calculate the M term on RHS of solution equation
    [~, M_linear, cw] = Mterm3D(mgrid, medium, medium_cond);
    
    % add non-rflecting layer
    M_linear = repmat(M_linear, 1, mgrid.num_x, mgrid.num_y) + M_gamma;
    
    % prepare for the calculation of wavevector along the propagation
    % direction
    Ktemp1_kx = zeros(1, mgrid.num_x);
    Ktemp1_ky = zeros(1,1,mgrid.num_y);

    Ktemp1_kx(1, :)  = squeeze(mgrid.kx);
    Ktemp1_ky(1,1,:) = squeeze(mgrid.ky);
    Ktemp = single(repmat((mgrid.w.'/medium.c0).^2,1,mgrid.num_x, mgrid.num_y) - ...
            repmat(Ktemp1_kx.^2, mgrid.num_t, 1, mgrid.num_y) -...
            repmat(Ktemp1_ky.^2, mgrid.num_t, mgrid.num_x, 1) - ...
            M_linear);

    K = single(sqrt(Ktemp));
    clear Ktemp
    if mod(mgrid.num_t,2) == 1
        K(mgrid.num_t/2+0.5:end,:,:) = -K(mgrid.num_t/2+0.5:end,:,:);
    else
        K(mgrid.num_t/2:end,:,:)     = -K(mgrid.num_t/2:end,:,:);
    end
    K(K==0)=eps;    
    
    f = excit_p;
    clear excit_p
    
    % Fourier transform of the excitation with regard of x, y and time
    excit_F = fftshift(fftn(f)); 
    h = waitbar(0,'Backward projection, please wait...');
    
    for I = mgrid.num_z:-1:1
        % used for filtering out evanescent wave
        Ktemp1 = repmat((mgrid.w.'./cw).^2,1,mgrid.num_x, mgrid.num_y)- ...
                 evanescent_fiter*repmat(Ktemp1_kx.^2, mgrid.num_t, 1, mgrid.num_y) -...
                 evanescent_fiter*repmat(Ktemp1_ky.^2, mgrid.num_t, mgrid.num_x, 1) - ...
                 evanescent_fiter*M_linear;
                  
        F5 = excit_F.*exp(1i*K*backward_dz); 
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
    
        excit_F = F5;
        f = real(ifftn(ifftshift(F5)));
        clear F5  
        
        if max(max(max(max(abs(f)))))==0
            disp ...
            ('The computation has lead to unphysical results. Please try reducing your step size')
            break
        end         
        
        if MDM_ani % show the animation
            p_ani(:,:, I) = squeeze(f(:,:,round(mgrid.num_y/2)));
            
            % in case a small computational time domain is applied
            it = round(I*mgrid.dz/medium.c0/mgrid.dt);
            
            imagesc (mgrid.z*1e3, mgrid.x*1e3, squeeze(p_ani(it,:,:)))
            axis image
            drawnow;  
            if MDM_mov % record the animation
                writeVideo(writerObj, getframe(gcf));
            end
        end
         
        if sum(sum(sum(sensor_mask(:,:,I))))~=0
            
            % calculate the starting position of recorded signal
            sensor_mask_I1 = sum(sum(sum(sensor_mask(:,:,1:I-1)))) + 1;
    
            % calculate the ending position of recorded signal
            sensor_mask_I2 = sensor_mask_I1 -1 +  ...
            sum(sum(sum(sensor_mask(:,:,I))));
    
            % save the time-domain signals
            p_time(:,sensor_mask_I1:sensor_mask_I2) = f(:,sensor_mask(:,:,I)~=0);
        end
            
        waitbar((mgrid.num_z - I+1)/mgrid.num_z)
    end
        
    close(h)

%%
elseif medium_cond == 2
    
    disp ('Homogeneous medium')
    
    % calculate the M term on RHS of solution equation
    [M_nonlinear, M_linear, cw] = Mterm3D(mgrid, medium, medium_cond);       
    
    Ktemp1_kx = zeros(1, mgrid.num_x);
    Ktemp1_ky = zeros(1,1,mgrid.num_y);

    Ktemp1_kx(1, :)  = single(squeeze(mgrid.kx));
    Ktemp1_ky(1,1,:) = single(squeeze(mgrid.ky));

    Ktemp = single(repmat((mgrid.w.'.^2),1,mgrid.num_x, mgrid.num_y)/medium.c0^2 - ...
            repmat(Ktemp1_kx.^2, mgrid.num_t, 1, mgrid.num_y) -...
            repmat(Ktemp1_ky.^2, mgrid.num_t, mgrid.num_x, 1) - ...
            repmat(M_linear, 1, mgrid.num_x, mgrid.num_y));

    K = single(sqrt(Ktemp));
    clear Ktemp
    if mod(mgrid.num_t,2) == 1
       K(mgrid.num_t/2+0.5:end,:,:) = -K(mgrid.num_t/2+0.5:end,:,:);
    else
       K(mgrid.num_t/2:end,:,:)     = -K(mgrid.num_t/2:end,:,:);
    end    
    
    K(K==0)=eps;
    
    % used for filtering out evanescent wave
    Ktemp1 = single(repmat(((mgrid.w'./cw).^2),1,mgrid.num_x, mgrid.num_y) - ...
             evanescent_fiter*repmat(Ktemp1_kx.^2, mgrid.num_t, 1, mgrid.num_y) -...
             evanescent_fiter*repmat(Ktemp1_ky.^2, mgrid.num_t, mgrid.num_x, 1) - ...
             evanescent_fiter*repmat(M_linear, 1, mgrid.num_x, mgrid.num_y));   
    clear M_linear    Ktemp1_kx    Ktemp1_ky   cw    
    
    f = excit_p;
    clear excit_p
    
    % Fourier transform of the excitation with regard of x, y and time
    excit_F = fftshift(fftn(f)); 

    h = waitbar(0,'Backward projection, please wait...');
    
    for I = mgrid.num_z:-1:1
        
        % Simpson 
        % 1   
        M  = repmat(M_nonlinear, 1, mgrid.num_x, mgrid.num_y).*...
             fftshift(fftn(f.*f));
        F1 = single(excit_F.*exp(0.5*1i*K*backward_dz) +...
             0.5*backward_dz*M.*exp(0.5*1i*K*backward_dz)./(2i.*K));  % P01(z+backward_dz/2)
        F1(isnan(F1)) = 0;
        F1(real(Ktemp1)<=0) = 0; 
    
        % 2
        clear f
        f = real(ifftn(ifftshift(F1)));  % p0(z+backward_dz/2) 
        clear F1
        M1 = repmat(M_nonlinear, 1, mgrid.num_x, mgrid.num_y).*...
             fftshift(fftn(f.*f));
        clear f
        F2 = excit_F.*exp(0.5*1i*K*backward_dz) +...
            0.25*backward_dz*exp(0.5*1i*K*backward_dz)./...
            (2i.*K).*(M + M1./exp(0.5*1i*K*backward_dz));% P12(z+backward_dz/2)
        clear M1
        F2(isnan(F2)) = 0;
        F2(real(Ktemp1)<=0) = 0;   
  
        % 3   
        f   = real(ifftn(ifftshift(F2))); % p1(z+backward_dz/2)   
        M2 = repmat(M_nonlinear, 1, mgrid.num_x, mgrid.num_y).*...
             fftshift(fftn(f.*f));
        clear f
        F3 = F2.*exp(0.5*1i*K*backward_dz) +...
             0.5*backward_dz*M2.*exp(0.5*1i*K*backward_dz)./(2i.*K);  % P03(z+backward_dz)
        F3(isnan(F3)) = 0;
        F3(real(Ktemp1)<=0) = 0;  
    
        % 4   
        f   = real(ifftn(ifftshift(F3))); % p0(z+backward_dz)  
        clear F3
        M3 = repmat(M_nonlinear, 1, mgrid.num_x, mgrid.num_y).*...
             fftshift(fftn(f.*f));
        clear f
        F4 = F2.*exp(0.5*1i*K*backward_dz) + ...
             0.25*backward_dz.*(M2 + M3./exp(0.5*1i*K*backward_dz)).*...
             exp(0.5*1i*K*backward_dz)./(2i.*K); % P14(z+backward_dz)
        clear M3   F2   
        F4(isnan(F4)) = 0;
        F4(real(Ktemp1)<=0) = 0;    
      
        % 5   
        f   = real(ifftn(ifftshift(F4))); % p0(z+backward_dz) 
        clear F4
        M4 = repmat(M_nonlinear, 1, mgrid.num_x, mgrid.num_y).*...
             fftshift(fftn(f.*f));
        clear f
        F5 = excit_F.*exp(1i*K*backward_dz) + ...
             backward_dz/6.0.*(M + 4*M2./exp(0.5*1i*K*backward_dz) +...
             M4./exp(1i*K*backward_dz)).*...
             exp(1i*K*backward_dz)./(2i.*K); % P25(z+backward_dz)
       
        clear M   M2   M4 
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
    
        excit_F = F5;
        f = real(ifftn(ifftshift(F5))); % p0(z+backward_dz) 
        clear F5

        if max(max(max(max(abs(f)))))==0
            disp ('The computation has lead to unphysical results. Please try reducing your step size.')
            break
        end 
         if MDM_ani % show the animation
             p_ani(:,:, I) = squeeze(f(:,:,round(mgrid.num_y/2)));
             
             % in case a small computational time domain is applied
             it = round(I*mgrid.dz/medium.c0/mgrid.dt);
             
             imagesc (mgrid.z*1e3, mgrid.x*1e3, squeeze(p_ani(it,:,:)))
             axis image
             drawnow;  
             if MDM_mov % record the animation
                 writeVideo(writerObj, getframe(gcf));
             end
         end  
         
        if sum(sum(sum(sensor_mask(:,:,I))))~=0
           % calculate the starting position of recorded signal
           sensor_mask_I1 = sum(sum(sum(sensor_mask(:,:,1:I-1)))) + 1;
    
           % calculate the ending position of recorded signal
           sensor_mask_I2 = sensor_mask_I1 -1 +  ...
           sum(sum(sum(sensor_mask(:,:,I))));
   
           % save the time-domain signals
           p_time(:,sensor_mask_I1:sensor_mask_I2) = f(:,sensor_mask(:,:,I)~=0);
        end
        
        waitbar((mgrid.num_z - I+1)/mgrid.num_z)
    end
    close(h)
%%    

elseif medium_cond == 3
    
    disp ('Inhomogeneous only in nonlinearity coefficient')
        
    % calculate the M term on RHS of solution equation
    [M_nonlinear, M_linear, cw] = Mterm3D(mgrid, medium, medium_cond);  
    
    % add absorption layer
    M_linear = repmat(M_linear, 1, mgrid.num_x, mgrid.num_y) + M_gamma;
    
    % prepare for the calculation of evanescent wave
    Ktemp1_kx = zeros(1, mgrid.num_x);
    Ktemp1_ky = zeros(1,1,mgrid.num_y);

    Ktemp1_kx(1, :)  = squeeze(mgrid.kx);
    Ktemp1_ky(1,1,:) = squeeze(mgrid.ky);

    Ktemp = repmat((mgrid.w.'./medium.c0).^2,1,mgrid.num_x, mgrid.num_y) - ...
            repmat(Ktemp1_kx.^2, mgrid.num_t, 1, mgrid.num_y) -...
            repmat(Ktemp1_ky.^2, mgrid.num_t, mgrid.num_x, 1) - ...
            M_linear;    
        
    K = sqrt(Ktemp);
    clear Ktemp
    if mod(mgrid.num_t,2) == 1
        K(mgrid.num_t/2+0.5:end,:,:) = -K(mgrid.num_t/2+0.5:end,:,:);
    else
        K(mgrid.num_t/2:end,:,:)     = -K(mgrid.num_t/2:end,:,:);
    end
    K(K==0)=eps;    
 
    % exponential term for the forward propagation
    expn2 = exp(0.5*1i*K*mgrid.dz);    
    
    Ktemp1 = repmat(((mgrid.w'./cw).^2),1,mgrid.num_x, mgrid.num_y)-...
             evanescent_fiter*repmat(Ktemp1_kx.^2, mgrid.num_t, 1, mgrid.num_y) -...
             evanescent_fiter*repmat(Ktemp1_ky.^2, mgrid.num_t, mgrid.num_x, 1) -...
             evanescent_fiter*M_linear;  
    clear  Ktemp1_c   Ktemp1_kx  Ktemp1_ky  M_linear
    
    % normalized wave field
    f = excit_p;
        
    % Fourier transform of the excitation with regard of x, y and time
    excit_F = fftshift(fftn(f)); 
    
    h = waitbar(0,'Backward projection, please wait...');

    for I = mgrid.num_z:-1:1
        
        % Simpson 
        % 1   
        F2t = fftshift(fft(f.*f,[],1),1);
        % M(f(z))
        M = fftshift(fft(fftshift(fft(M_nonlinear(:,:,:,I+1).*F2t,[],2),2),[],3),3);
        clear   Ft  F2t
        F1 = excit_F.*expn2 + 0.5*backward_dz*M.*expn2./(2i.*K);  % P01(z+backward_dz/2)
        F1(isnan(F1)) = 0;
        F1(real(Ktemp1)<=0) = 0; 
    
        % 2
        f = real(ifftn(ifftshift(F1)));  % p0(z+backward_dz/2) 
        clear F1
        F2t = fftshift(fft(f.*f,[],1),1);  % fft(f^2)   
        % M(f0(z+backward_dz/2))  M1
        M1 = fftshift(fft(fftshift(fft(M_nonlinear(:,:,:,I+1).*F2t,[],2),2),[],3),3);
        clear   Ft  F2t
        F2 = excit_F.*expn2 + 0.25*backward_dz*expn2./(2i.*K).*(M + M1./expn2);% P12(z+backward_dz/2)
        clear M1
        F2(isnan(F2)) = 0;
        F2(real(Ktemp1)<=0) = 0;   
        
        % 3   
        f   = real(ifftn(ifftshift(F2))); % p1(z+backward_dz/2)   
        F2t = fftshift(fft(f.*f,[],1),1); % fft(f1^2)   
        % M(f1(z+backward_dz/2)) M2
        M2 = fftshift(fft(fftshift(fft(M_nonlinear(:,:,:,I+1).*F2t,[],2),2),[],3),3);
        clear   Ft  F2t
        F3 = F2.*expn2 + 0.5*backward_dz*M2.*expn2./(2i.*K);  % P03(z+backward_dz)
        F3(isnan(F3)) = 0;
        F3(real(Ktemp1)<=0) = 0;  
    
        % 4     
        f   = real(ifftn(ifftshift(F3))); % p0(z+backward_dz)   
        F2t = fftshift(fft(f.*f,[],1),1); % fft(f2^2)
        % M(f0(z+backward_dz))  M3
        M3 = fftshift(fft(fftshift(fft(M_nonlinear(:,:,:,I).*F2t,[],2),2),[],3),3);
        F4 = F2.*expn2 + 0.25*backward_dz.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z+backward_dz)
        clear M3   Ft  F2t
        F4(isnan(F4)) = 0;
        F4(real(Ktemp1)<=0) = 0;    
      
        % 5   
        f   = real(ifftn(ifftshift(F4))); % p0(z+backward_dz) 
        clear F4
        F2t = fftshift(fft(f.*f,[],1),1); % fft(f2^2)
        % M(f1(z+backward_dz)) M4
        M4 = fftshift(fft(fftshift(fft(M_nonlinear(:,:,:,I).*F2t,[],2),2),[],3),3);
        F5 = excit_F.*exp(1i*K*backward_dz) + ...
             backward_dz/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*backward_dz)).*...
             exp(1i*K*backward_dz)./(2i.*K); % P25(z+backward_dz)
     
        clear M   M2   M4  Ft  F2t
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
    
        excit_F = F5;
        f = real(ifftn(ifftshift(F5))); % p0(z+backward_dz) 
        clear F5
        
        if max(max(max(max(abs(f)))))==0
            disp ('The computation has lead to unphysical results. Please try reducing your step size.')
            break
        end   
         if MDM_ani % show the animation
             p_ani(:,:, I) = squeeze(f(:,:,round(mgrid.num_y/2)));
             
             % in case a small computational time domain is applied
             it = round(I*mgrid.dz/medium.c0/mgrid.dt);
             
             imagesc (mgrid.z*1e3, mgrid.x*1e3, squeeze(p_ani(it,:,:)))
             axis image
             drawnow;  
             if MDM_mov % record the animation
                 writeVideo(writerObj, getframe(gcf));
             end
         end  

        if sum(sum(sum(sensor_mask(:,:,I))))~=0
           % calculate the starting position of recorded signal
           sensor_mask_I1 = sum(sum(sum(sensor_mask(:,:,1:I-1)))) + 1;
    
           % calculate the ending position of recorded signal
           sensor_mask_I2 = sensor_mask_I1 -1 +  ...
           sum(sum(sum(sensor_mask(:,:,I))));
    
           % save the time-domain signals
           p_time(:,sensor_mask_I1:sensor_mask_I2) = f(:,sensor_mask(:,:,I)~=0);
        end        
        
        waitbar((mgrid.num_z - I+1)/mgrid.num_z)
    end
    close(h)
        

elseif medium_cond == 4 && MDMC == 0
    
    disp ('linear media')
    % prepare for the calculation of evanescent wave
    Ktemp1_kx = zeros(1, mgrid.num_x);
    Ktemp1_ky = zeros(1,1,mgrid.num_y);

    Ktemp1_kx(1, :)  = squeeze(mgrid.kx);
    Ktemp1_ky(1,1,:) = squeeze(mgrid.ky);

    Ktemp = repmat((mgrid.w.'./medium.c0).^2,1,mgrid.num_x, mgrid.num_y) - ...
            repmat(Ktemp1_kx.^2, mgrid.num_t, 1, mgrid.num_y) -...
            repmat(Ktemp1_ky.^2, mgrid.num_t, mgrid.num_x, 1);
     
    K = sqrt(Ktemp);
    clear Ktemp
    if mod(mgrid.num_t,2) == 1
        K(mgrid.num_t/2+0.5:end,:,:) = -K(mgrid.num_t/2+0.5:end,:,:);
    else
        K(mgrid.num_t/2:end,:,:)     = -K(mgrid.num_t/2:end,:,:);
    end
    K(K==0)=eps;        
        
    rho_rho = medium.rho;
    if length(medium.rho) == 1
        rho_rho = medium.rho*ones(mgrid.num_x, mgrid.num_y, mgrid.num_z+1);
    end       
    
    % exponential term for the forward propagation
    expn2 = exp(0.5*1i*K*mgrid.dz);         
    
    rho2 = zeros(1, mgrid.num_x, mgrid.num_y);   
    rho2(1,:,:) = squeeze(rho_rho(:,:,end));
        
    % normalized wave field
    f = excit_p./sqrt(repmat(rho2, mgrid.num_t, 1,1));

    % Fourier transform of the excitation with regard of x, y and time
    excit_F = fftshift(fftn(f)); 
        
    % calculate the M term on RHS of solution equation
    [~, M_linear, cw] = Mterm3D(mgrid, medium, medium_cond); 
    
    M_nonlinear = 0;
    % add absorption layer
    M_linear = M_linear + repmat(M_gamma, 1,1,1,mgrid.num_z+1);   
    
    h = waitbar(0,'Backward projection, please wait...');
        
    for I = mgrid.num_z:-1:1
        
        % used for filtering out evanescent wave
        Ktemp1 = (repmat((mgrid.w.'),1,mgrid.num_x, mgrid.num_y)./...
                 cw(:,:,:,I)).^2 - ...
                 evanescent_fiter*repmat(Ktemp1_kx.^2, mgrid.num_t, 1, mgrid.num_y) -...
                 evanescent_fiter*repmat(Ktemp1_ky.^2, mgrid.num_t, mgrid.num_x, 1);
                 
        % Simpson 
        % 1   
        Ft  = fftshift(fft(f,[],1),1);
        % M(f(z))
        M = fftshift(fft(fftshift(fft((M_linear(:,:,:,I+1)).*Ft,[],2),2),[],3),3);
        clear   Ft  F2t
        F1 = excit_F.*expn2 + 0.5*backward_dz*M.*expn2./(2i.*K);  % P01(z+backward_dz/2)
        F1(isnan(F1)) = 0;
        F1(real(Ktemp1)<=0) = 0; 
   
        % 2
        f = real(ifftn(ifftshift(F1)));  % p0(z+backward_dz/2) 
        clear F1
        Ft  = fftshift(fft(f,[],1),1);     % fft(f) 
        % M(f0(z+backward_dz/2))  M1
        M1 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I+1)).*Ft,[],2),2),[],3),3);
        clear   Ft  F2t
        F2 = excit_F.*expn2 + 0.25*backward_dz*expn2./(2i.*K).*(M + M1./expn2);% P12(z+backward_dz/2)
        clear M1
        F2(isnan(F2)) = 0;
        F2(real(Ktemp1)<=0) = 0;   
  
        % 3   
        f   = real(ifftn(ifftshift(F2))); % p1(z+backward_dz/2)   
        Ft  = fftshift(fft(f,[],1),1);    % fft(f1) 
        % M(f1(z+backward_dz/2)) M2
        M2 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I+1)).*Ft,[],2),2),[],3),3);
        clear   Ft  F2t
        F3 = F2.*expn2 + 0.5*backward_dz*M2.*expn2./(2i.*K);  % P03(z+backward_dz)
        F3(isnan(F3)) = 0;
        F3(real(Ktemp1)<=0) = 0;  
    
        % 4     
        f   = real(ifftn(ifftshift(F3))); % p0(z+backward_dz)   
        Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
        % M(f0(z+backward_dz))  M3
        M3 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I)).*Ft,[],2),2),[],3),3);
        F4 = F2.*expn2 + 0.25*backward_dz.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z+backward_dz)
        clear M3   Ft  F2t
        F4(isnan(F4)) = 0;
        F4(real(Ktemp1)<=0) = 0;    
      
        % 5   
        f   = real(ifftn(ifftshift(F4))); % p0(z+backward_dz) 
        clear F4
        Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
        % M(f1(z+backward_dz)) M4
        M4 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I)).*Ft,[],2),2),[],3),3);
        F5 = excit_F.*exp(1i*K*backward_dz) + ...
             backward_dz/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*backward_dz)).*...
             exp(1i*K*backward_dz)./(2i.*K); % P25(z+backward_dz)
     
        clear M   M2   M4  Ft  F2t
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
         
        excit_F = F5;
        f = real(ifftn(ifftshift(F5))); % p0(z+backward_dz) 
        clear F5

        if max(max(max(max(abs(f)))))==0
            disp ('The computation has lead to unphysical results. Please try reducing your step size.')
            break
        end         
        
%% reflection
        rho2(1,:,:) = rho_rho(:,:,I);
        if reflection_order~=0
            Ktemp2 = repmat((mgrid.w.'.^2),1,mgrid.num_x, mgrid.num_y)./...
                     squeeze(cw(:,:,:,I+1)).^2 - ...
                     repmat(Ktemp1_kx.^2, mgrid.num_t, 1, mgrid.num_y) -...
                     repmat(Ktemp1_ky.^2, mgrid.num_t, mgrid.num_x, 1);                          
            K2 = sqrt(Ktemp2);
            if mod(mgrid.num_t,2) == 1
                K2(mgrid.num_t/2+0.5:end,:,:) = -K2(mgrid.num_t/2+0.5:end,:,:);
            else
                K2(mgrid.num_t/2:end,:,:)    = -K2(mgrid.num_t/2:end,:,:);
            end
            K2(K2==0) = eps;
            clear Ktemp2            
           
            rho1(1,:,:) = rho_rho(:,:,I+1);
            rho11 = repmat(rho1, mgrid.num_t, 1,1); 
            rho22 = repmat(rho2, mgrid.num_t, 1,1);            
            T_rhoc = 2.*rho11.*cw(:,:,:,I+1)./...
                    (rho11.*cw(:,:,:,I+1) + rho22.*cw(:,:,:,I));          
            clear rho11 rho22 rho1   
            P_in = excit_F.*exp(1i.*K2*mgrid.dz);
            p_ref(:,:,:, I) = real(ifftn(ifftshift(P_in))).*(T_rhoc-1);
            clear P_in  T_rhoc  K2 
        end        

        % recover the normalized wave field
        p_I = sqrt(repmat(rho2, mgrid.num_t, 1,1)).*f;    
        
         if MDM_ani % show the animation
             p_ani(:,:,I) = squeeze(p_I(:,:,round(mgrid.num_y/2)));
             
             % in case a small computational time domain is applied
             it = round(I*mgrid.dz/medium.c0/mgrid.dt);
             
             imagesc (mgrid.z*1e3, mgrid.x*1e3, squeeze(p_ani(it,:,:)))
             axis image
             drawnow;  
             if MDM_mov % record the animation
                 writeVideo(writerObj, getframe(gcf));
             end
         end    
         
        if sum(sum(sum(sensor_mask(:,:,I))))~=0
            % calculate the starting position of recorded signal
            sensor_mask_I1 = sum(sum(sum(sensor_mask(:,:,1:I-1)))) + 1;
            % calculate the ending position of recorded signal
            sensor_mask_I2 = sensor_mask_I1 -1 +  ...
            sum(sum(sum(sensor_mask(:,:,I))));
    
            % save the time-domain signals
            p_time(:,sensor_mask_I1:sensor_mask_I2) = p_I(:,sensor_mask(:,:,I)~=0);
        end
        clear p_I
        
        waitbar((mgrid.num_z - I+1)/mgrid.num_z)
    end
    close(h)

        
elseif medium_cond ==4 && MDMC == 1 % add correction  
    disp ('Linear media with transmission corrections')
    
    % prepare for the calculation of wavevector along the propagation
    % direction
    Ktemp1_kx = zeros(1, mgrid.num_x);
    Ktemp1_ky = zeros(1,1,mgrid.num_y);

    Ktemp1_kx(1, :)  = squeeze(mgrid.kx);
    Ktemp1_ky(1,1,:) = squeeze(mgrid.ky);
    
    Ktemp = repmat((mgrid.w.'./medium.c0).^2,1,mgrid.num_x, mgrid.num_y) - ...
            repmat(Ktemp1_kx.^2, mgrid.num_t, 1, mgrid.num_y) -...
            repmat(Ktemp1_ky.^2, mgrid.num_t, mgrid.num_x, 1);
     
    K = sqrt(Ktemp);
    clear Ktemp
    if mod(mgrid.num_t,2) == 1
        K(mgrid.num_t/2+0.5:end,:,:) = -K(mgrid.num_t/2+0.5:end,:,:);
    else
        K(mgrid.num_t/2:end,:,:)     = -K(mgrid.num_t/2:end,:,:);
    end
    K(K==0)=eps;        
        
    rho_rho = medium.rho;
    if length(medium.rho) == 1
        rho_rho = medium.rho*ones(mgrid.num_x, mgrid.num_y, mgrid.num_z+1);
    end      
    
    % exponential term for the forward propagation
    expn2 = exp(0.5*1i*K*backward_dz);         
    
    rho1 = zeros(1, mgrid.num_x, mgrid.num_y);
    rho2 = rho1;    

    % normalized wave field
    f = excit_p;
    clear excit_p
    
    % Fourier transform of the excitation with regard of x, y and time
    excit_F = fftshift(fftn(f)); 
        
    % calculate the M term on RHS of solution equation
    [~, M_linear, cw] = Mterm3D_MMDM(mgrid, medium, medium_cond);   
    
    M_linear = M_linear + repmat(M_gamma, 1,1,1,mgrid.num_z+1);  
    M_nonlinear = 0;
    
    h = waitbar(0,'Backward projection, please wait...');
        
    for I = mgrid.num_z:-1:1
 
        Ktemp2 = (repmat((mgrid.w.'),1,mgrid.num_x, mgrid.num_y)./...
                 cw(:,:,:,I)).^2 - ...
                 repmat(Ktemp1_kx.^2, mgrid.num_t, 1, mgrid.num_y) -...
                 repmat(Ktemp1_ky.^2, mgrid.num_t, mgrid.num_x, 1);                          
        K2 = sqrt(Ktemp2);
        if mod(mgrid.num_t,2) == 1
            K2(mgrid.num_t/2+0.5:end,:,:) = -K2(mgrid.num_t/2+0.5:end,:,:);
        else
            K2(mgrid.num_t/2:end,:,:)    = -K2(mgrid.num_t/2:end,:,:);
        end
        K2(K2==0) = eps;
        clear Ktemp2  
        
        % used for filtering out evanescent wave
        Ktemp1 = (repmat((mgrid.w.'),1,mgrid.num_x, mgrid.num_y)./...
                 cw(:,:,:,I)).^2 - ...
                 evanescent_fiter*repmat(Ktemp1_kx.^2, mgrid.num_t, 1, mgrid.num_y) -...
                 evanescent_fiter*repmat(Ktemp1_ky.^2, mgrid.num_t, mgrid.num_x, 1);
       
        % Simpson 
        % 1   
        Ft  = fftshift(fft(f,[],1),1);
        % M(f(z))
        M = fftshift(fft(fftshift(fft((M_linear(:,:,:,I+1)).*Ft,[],2),2),[],3),3);
        clear   Ft  F2t
        F1 = excit_F.*expn2 + 0.5*backward_dz*M.*expn2./(2i.*K);  % P01(z+backward_dz/2)
        F1(isnan(F1)) = 0;
        F1(real(Ktemp1)<=0) = 0; 
   
        % 2
        f = real(ifftn(ifftshift(F1)));  % p0(z+backward_dz/2) 
        clear F1
        Ft  = fftshift(fft(f,[],1),1);     % fft(f) 
        % M(f0(z+backward_dz/2))  M1
        M1 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I+1)).*Ft,[],2),2),[],3),3);
        clear   Ft  F2t
        F2 = excit_F.*expn2 + 0.25*backward_dz*expn2./(2i.*K).*(M + M1./expn2);% P12(z+backward_dz/2)
        clear M1
        F2(isnan(F2)) = 0;
        F2(real(Ktemp1)<=0) = 0;   
  
        % 3   
        f   = real(ifftn(ifftshift(F2))); % p1(z+backward_dz/2)   
        Ft  = fftshift(fft(f,[],1),1);    % fft(f1) 
        % M(f1(z+backward_dz/2)) M2
        M2 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I+1)).*Ft,[],2),2),[],3),3);
        clear   Ft  F2t
        F3 = F2.*expn2 + 0.5*backward_dz*M2.*expn2./(2i.*K);  % P03(z+backward_dz)
        F3(isnan(F3)) = 0;
        F3(real(Ktemp1)<=0) = 0;  
    
        % 4     
        f   = real(ifftn(ifftshift(F3))); % p0(z+backward_dz)   
        Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
        % M(f0(z+backward_dz))  M3
        M3 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I)).*Ft,[],2),2),[],3),3);
        F4 = F2.*expn2 + 0.25*backward_dz.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z+backward_dz)
        clear M3   Ft  F2t
        F4(isnan(F4)) = 0;
        F4(real(Ktemp1)<=0) = 0;    
      
        % 5   
        f   = real(ifftn(ifftshift(F4))); % p0(z+backward_dz) 
        clear F4
        Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
        % M(f1(z+backward_dz)) M4
        M4 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I)).*Ft,[],2),2),[],3),3);
        F5 = excit_F.*exp(1i*K*backward_dz) + ...
             backward_dz/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*backward_dz)).*...
             exp(1i*K*backward_dz)./(2i.*K); % P25(z+backward_dz)
     
        clear M   M2   M4  Ft  F2t
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
  
%% add phase correction
        wxy2   = repmat(mgrid.w.', 1,mgrid.num_x, mgrid.num_y).*...
                 repmat(mgrid.w.', 1,mgrid.num_x, mgrid.num_y);
        e1     = exp(1i*(K - K2)*backward_dz);
        k_dif2 = (wxy2./(medium.c0*medium.c0) - wxy2./(cw(:,:,:,I).*cw(:,:,:,I)))*backward_dz;        
        k_dif  = (wxy2./(medium.c0*medium.c0) - wxy2./(cw(:,:,:,I+1).*cw(:,:,:,I+1)))*backward_dz;         
        
        correction = excit_F.*exp(1i*K2*backward_dz).*...
                    (1-e1.*(1+k_dif./4i./K)./(1-k_dif2./4i./K));
        clear k_dif  k_dif2  el  wxy2        
        
        F5 = F5 + correction;  
        clear correction
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
        f = real(ifftn(ifftshift(F5)));
        
        rho1(1,:,:) = rho_rho(:,:,I+1);
        rho2(1,:,:) = rho_rho(:,:, I);        
        rho11 = repmat(rho1, mgrid.num_t, 1,1); 
        rho22 = repmat(rho2, mgrid.num_t, 1,1);            
        T_rhoc = 2.*rho11.*cw(:,:,:,I+1)./...
                (rho11.*cw(:,:,:,I+1) + rho22.*cw(:,:,:,I));          
        clear rho11 rho22 rho1         
        f = f./T_rhoc;
        F5 = fftshift(fftn(f));
        
%%        
        excit_F = F5;
        f = real(ifftn(ifftshift(F5))); % p0(z+backward_dz) 
        clear F5
        
        if max(max(max(max(abs(f)))))==0
            disp ('The computation has lead to unphysical results. Please try reducing your step size.')
            break
        end         
        
        if reflection_order~=0
            Ktemp3 = (repmat(mgrid.w.',1,mgrid.num_x, mgrid.num_y)./...
                     cw(:,:,:,I+1)).^2 - ...
                     repmat(Ktemp1_kx.^2, mgrid.num_t, 1, mgrid.num_y) -...
                     repmat(Ktemp1_ky.^2, mgrid.num_t, mgrid.num_x, 1);                          
            K3 = sqrt(Ktemp3);
            if mod(mgrid.num_t,2) == 1
                K3(mgrid.num_t/2+0.5:end,:,:) = -K3(mgrid.num_t/2+0.5:end,:,:);
            else
            K3(mgrid.num_t/2:end,:,:)    = -K3(mgrid.num_t/2:end,:,:);
            end
            K3(K3==0) = eps;
            clear Ktemp3              
            P_in = excit_F.*exp(1i.*K3*mgrid.dz);
            p_ref(:,:,:,I) = real(ifftn(ifftshift(P_in))).*(T_rhoc-1);
        end        
        
         if MDM_ani % show the animation
             p_ani(:,:,I) = squeeze(f(:,:,round(mgrid.num_y/2)));
             
             it = round(I*mgrid.dz/medium.c0/mgrid.dt);
             
             imagesc (mgrid.z*1e3, mgrid.x*1e3, squeeze(p_ani(it,:,:)))
             axis image
             drawnow;  
             if MDM_mov % record the animation
                 writeVideo(writerObj, getframe(gcf));
             end
         end  
         
        if sum(sum(sum(sensor_mask(:,:,I))))~=0
           % calculate the starting position of recorded signal
           sensor_mask_I1 = sum(sum(sum(sensor_mask(:,:,1:I-1)))) + 1;
    
           % calculate the ending position of recorded signal
           sensor_mask_I2 = sensor_mask_I1 -1 +  ...
           sum(sum(sum(sensor_mask(:,:,I))));
    
           % save the time-domain signals
           p_time(:,sensor_mask_I1:sensor_mask_I2) = f(:,sensor_mask(:,:,I)~=0);
        end
        clear p_I
        
        waitbar((mgrid.num_z - I+1)/mgrid.num_z)
    end
    close(h)    
        
elseif medium_cond == 5 && MDMC ==0
    
%     disp ('Nonlinear media without tranmission corrections')
    
    % calculate the M term on RHS of solution equation
    [M_nonlinear, M_linear, cw] = Mterm3D(mgrid, medium, medium_cond);    
    
    % add absorption layers
    M_linear = M_linear + repmat(M_gamma, 1,1,1,mgrid.num_z+1);  
    
    % prepare for the calculation of evanescent wave
    Ktemp1_kx = zeros(1, mgrid.num_x);
    Ktemp1_ky = zeros(1,1,mgrid.num_y);

    Ktemp1_kx(1, :)  = squeeze(mgrid.kx);
    Ktemp1_ky(1,1,:) = squeeze(mgrid.ky);

    Ktemp = repmat((mgrid.w.'./medium.c0).^2,1,mgrid.num_x, mgrid.num_y) - ...
            repmat(Ktemp1_kx.^2, mgrid.num_t, 1, mgrid.num_y) -...
            repmat(Ktemp1_ky.^2, mgrid.num_t, mgrid.num_x, 1);    
    
    K = sqrt(Ktemp);
    clear Ktemp
    if mod(mgrid.num_t,2) == 1
        K(mgrid.num_t/2+0.5:end,:,:) = -K(mgrid.num_t/2+0.5:end,:,:);
    else
        K(mgrid.num_t/2:end,:,:)     = -K(mgrid.num_t/2:end,:,:);
    end
    K(K==0)=eps;    
    
    rho_rho = medium.rho;
    if length(medium.rho) == 1
        rho_rho = medium.rho*ones(mgrid.num_x, mgrid.num_y, mgrid.num_z+1);
    end      

    rho2 = zeros(1, mgrid.num_x, mgrid.num_y);
    rho2(1,:,:) = squeeze(rho_rho(:,:,end));
    
    % normalized wave field
    f = excit_p./sqrt(repmat(rho2, mgrid.num_t, 1,1));
    
    % Fourier transform of the excitation with regard of x, y and time
    excit_F = fftshift(fftn(f)); 

    % exponential term for the forward propagation
    expn2 = exp(0.5*1i*K*backward_dz); 

    h = waitbar(0,'Backward projection, please wait...');

    for I = mgrid.num_z:-1:1
         
         % used for filtering out evanescent wave
         Ktemp1 = (repmat((mgrid.w.'),1,mgrid.num_x, mgrid.num_y)./...
                  cw(:,:,:,I)).^2 - ...
                  evanescent_fiter*repmat(Ktemp1_kx.^2, mgrid.num_t, 1, mgrid.num_y) -...
                  evanescent_fiter*repmat(Ktemp1_ky.^2, mgrid.num_t, mgrid.num_x, 1);
         
         % Simpson 

         % 1   
         Ft  = fftshift(fft(f,[],1),1);
         F2t = fftshift(fft(f.*f,[],1),1);
         clear f
         % M(f(z))
         M = fftshift(fft(fftshift(fft((M_linear(:,:,:,I+1)).*Ft,[],2),2),[],3),3);
         M = M + fftshift(fft(fftshift(fft(M_nonlinear(:,:,:,I+1).*F2t,[],2),2),[],3),3);
         clear   Ft  F2t
         F1 = excit_F.*expn2 + 0.5*backward_dz*M.*expn2./(2i.*K);  % P01(z+backward_dz/2)
         F1(isnan(F1)) = 0;
         F1(real(Ktemp1)<=0) = 0; 
    
         % 2
         f = real(ifftn(ifftshift(F1)));  % p0(z+backward_dz/2) 
         clear F1
         Ft  = fftshift(fft(f,[],1),1);     % fft(f) 
         F2t = fftshift(fft(f.*f,[],1),1);  % fft(f^2)   
         clear f
         % M(f0(z+backward_dz/2))  M1
         M1 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I+1)).*Ft,[],2),2),[],3),3);
         M1 = M1 + fftshift(fft(fftshift(fft(M_nonlinear(:,:,:,I+1).*F2t,[],2),2),[],3),3);
         clear   Ft  F2t
         F2 = excit_F.*expn2 + 0.25*backward_dz*expn2./(2i.*K).*(M + M1./expn2);% P12(z+backward_dz/2)
         clear M1
         F2(isnan(F2)) = 0;
         F2(real(Ktemp1)<=0) = 0;   
  
         % 3   
         f   = real(ifftn(ifftshift(F2))); % p1(z+backward_dz/2)   
         Ft  = fftshift(fft(f,[],1),1);    % fft(f1) 
         F2t = fftshift(fft(f.*f,[],1),1); % fft(f1^2)   
         clear f
         % M(f1(z+backward_dz/2)) M2
         M2 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I+1)).*Ft,[],2),2),[],3),3);
         M2 = M2 + fftshift(fft(fftshift(fft(M_nonlinear(:,:,:,I+1).*F2t,[],2),2),[],3),3);
         clear   Ft  F2t
         F3 = F2.*expn2 + 0.5*backward_dz*M2.*expn2./(2i.*K);  % P03(z+backward_dz)
         F3(isnan(F3)) = 0;
         F3(real(Ktemp1)<=0) = 0;  
    
         % 4     
         f   = real(ifftn(ifftshift(F3))); % p0(z+backward_dz)   
         Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
         F2t = fftshift(fft(f.*f,[],1),1); % fft(f2^2)
         clear f
         % M(f0(z+backward_dz))  M3
         M3 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I)).*Ft,[],2),2),[],3),3);
         M3 = M3 + fftshift(fft(fftshift(fft(M_nonlinear(:,:,:,I).*F2t,[],2),2),[],3),3);
         F4 = F2.*expn2 + 0.25*backward_dz.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z+backward_dz)
         clear M3   Ft  F2t
         F4(isnan(F4)) = 0;
         F4(real(Ktemp1)<=0) = 0;    
      
         % 5   
         f   = real(ifftn(ifftshift(F4))); % p0(z+backward_dz) 
         clear F4
         Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
         F2t = fftshift(fft(f.*f,[],1),1); % fft(f2^2)
         clear f
         % M(f1(z+backward_dz)) M4
         M4 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I)).*Ft,[],2),2),[],3),3);
         M4 = M4 + fftshift(fft(fftshift(fft(M_nonlinear(:,:,:,I).*F2t,[],2),2),[],3),3);
         F5 = excit_F.*exp(1i*K*backward_dz) + ...
              backward_dz/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*backward_dz)).*...
              exp(1i*K*backward_dz)./(2i.*K); % P25(z+backward_dz)
       
         clear M   M2   M4  Ft  F2t
         F5(isnan(F5)) = 0;
         F5(real(Ktemp1)<=0) = 0; 
 
 %%        
         excit_F = F5;
         f = real(ifftn(ifftshift(F5))); % p0(z+backward_dz) 
         clear F5
    
        if max(max(max(max(abs(f)))))==0
            disp ('The computation has lead to unphysical results. Please try reducing your step size.')
            break
        end 

        rho2(1,:,:) = squeeze(rho_rho(:,:, I));
        if reflection_order~=0
            Ktemp2 = (repmat((mgrid.w.'),1,mgrid.num_x, mgrid.num_y)./...
                     cw(:,:,:,I+1)).^2 - ...
                     repmat(Ktemp1_kx.^2, mgrid.num_t, 1, mgrid.num_y) -...
                     repmat(Ktemp1_ky.^2, mgrid.num_t, mgrid.num_x, 1);                          
            K2 = sqrt(Ktemp2);
            if mod(mgrid.num_t,2) == 1
                K2(mgrid.num_t/2+0.5:end,:,:) = -K2(mgrid.num_t/2+0.5:end,:,:);
            else
            K2(mgrid.num_t/2:end,:,:)    = -K2(mgrid.num_t/2:end,:,:);
            end
            K2(K2==0) = eps;
            clear Ktemp2            
 
            rho1(1,:,:) = rho_rho(:,:,I+1);
            rho11 = repmat(rho1, mgrid.num_t, 1,1); 
            rho22 = repmat(rho2, mgrid.num_t, 1,1);            
            T_rhoc = 2.*rho11.*cw(:,:,:,I+1)./...
                    (rho11.*cw(:,:,:,I+1) + rho22.*cw(:,:,:,I));          
            clear rho11 rho22 rho1                    
 
            P_in = excit_F.*exp(1i.*K2*mgrid.dz);
            p_ref(:,:,:,I) = real(ifftn(ifftshift(P_in))).*(T_rhoc-1);    
            clear T_rhoc  P_in  K2
        end   
        
         % recover the normalized wave field
         p_I = sqrt(repmat(rho2, mgrid.num_t, 1,1)).*f;    
         
         if MDM_ani % show the animation
             p_ani(:,:, I) = squeeze(p_I(:,:,round(mgrid.num_y/2)));
             
             % in case a small computational time domain is used
             it = round(I*mgrid.dz/medium.c0/mgrid.dt);
             
             imagesc (mgrid.z*1e3, mgrid.x*1e3, squeeze(p_ani(it,:,:)))
             axis image
             drawnow;  
             if MDM_mov % record the animation
                 writeVideo(writerObj, getframe(gcf));
             end
         end  
    
         if sum(sum(sum(sensor_mask(:,:,I))))~=0
             % calculate the starting position of recorded signal
             sensor_mask_I1 = sum(sum(sum(sensor_mask(:,:,1:I-1)))) + 1;
    
             % calculate the ending position of recorded signal
             sensor_mask_I2 = sensor_mask_I1 -1 +  ...
             sum(sum(sum(sensor_mask(:,:,I))));
    
             % save the time-domain signals
             p_time(:,sensor_mask_I1:sensor_mask_I2) = p_I(:,sensor_mask(:,:,I)~=0);
         end
         
         clear p_I
         
         waitbar((mgrid.num_z - I+1)/mgrid.num_z)
     end
    close(h)
    
%%

elseif medium_cond == 5 && MDMC == 1

%     disp ('Nonlinear media with transmission corrections')        
    % calculate the M term on RHS of solution equation
    [M_nonlinear, M_linear, cw] = Mterm3D_MMDM(mgrid, medium, medium_cond);            
     M_linear = M_linear + repmat(M_gamma, 1,1,1,mgrid.num_z+1);  
     
    % prepare for the calculation of evanescent wave
    Ktemp1_kx = zeros(1, mgrid.num_x);
    Ktemp1_ky = zeros(1,1,mgrid.num_y);

    Ktemp1_kx(1, :)  = squeeze(mgrid.kx);
    Ktemp1_ky(1,1,:) = squeeze(mgrid.ky);

    Ktemp = repmat((mgrid.w.'.^2),1,mgrid.num_x, mgrid.num_y)/medium.c0^2 - ...
            repmat(Ktemp1_kx.^2, mgrid.num_t, 1, mgrid.num_y) -...
            repmat(Ktemp1_ky.^2, mgrid.num_t, mgrid.num_x, 1);    
    
    K = sqrt(Ktemp);
    clear Ktemp
    if mod(mgrid.num_t,2) == 1
        K(mgrid.num_t/2+0.5:end,:,:) = -K(mgrid.num_t/2+0.5:end,:,:);
    else
        K(mgrid.num_t/2:end,:,:)     = -K(mgrid.num_t/2:end,:,:);
    end
    K(K==0)=eps;           

    rho_rho = medium.rho;
    if length(medium.rho) == 1
        rho_rho = medium.rho*ones(mgrid.num_x, mgrid.num_y, mgrid.num_z+1);
    end        
    
    rho1 = zeros(1, mgrid.num_x, mgrid.num_y);
    rho2 = zeros(1, mgrid.num_x, mgrid.num_y);   
    
    % normalized wave field
    f = excit_p;    
    clear  excit_p
    
    % Fourier transform of the excitation with regard of x, y and time
    excit_F = fftshift(fftn(f)); 

    % exponential term for the forward propagation
    expn2 = exp(0.5*1i*K*backward_dz);     
    
    h = waitbar(0,'Backward projection, please wait...');

    for I = mgrid.num_z:-1:1
         
        Ktemp2 = (repmat((mgrid.w.'),1,mgrid.num_x, mgrid.num_y)./...
                 cw(:,:,:,I)).^2 - ...
                 repmat(Ktemp1_kx.^2, mgrid.num_t, 1, mgrid.num_y) -...
                 repmat(Ktemp1_ky.^2, mgrid.num_t, mgrid.num_x, 1);                          
        K2 = sqrt(Ktemp2);
        if mod(mgrid.num_t,2) == 1
            K2(mgrid.num_t/2+0.5:end,:,:) = -K2(mgrid.num_t/2+0.5:end,:,:);
        else
            K2(mgrid.num_t/2:end,:,:) = -K2(mgrid.num_t/2:end,:,:);
        end
        K2(K2==0) = eps;
        clear Ktemp2                   
        
        % used for filtering out evanescent wave
        Ktemp1 = (repmat((mgrid.w.'),1,mgrid.num_x, mgrid.num_y)./...
                 cw(:,:,:,I)).^2 - ...
                 evanescent_fiter*repmat(Ktemp1_kx.^2, mgrid.num_t, 1, mgrid.num_y) -...
                 evanescent_fiter*repmat(Ktemp1_ky.^2, mgrid.num_t, mgrid.num_x, 1);

        % Simpson 
        % 1   
        Ft  = fftshift(fft(f,[],1),1);
        F2t = fftshift(fft(f.*f,[],1),1);
        clear f
        % M(f(z))
        M = fftshift(fft(fftshift(fft((M_linear(:,:,:,I+1)).*Ft,[],2),2),[],3),3);
        M = M + fftshift(fft(fftshift(fft(M_nonlinear(:,:,:,I+1).*F2t,[],2),2),[],3),3);
        clear   Ft  F2t
        F1 = excit_F.*expn2 + 0.5*backward_dz*M.*expn2./(2i.*K);  % P01(z+backward_dz/2)
        F1(isnan(F1)) = 0;
        F1(real(Ktemp1)<=0) = 0; 
    
        % 2
        f = real(ifftn(ifftshift(F1)));  % p0(z+backward_dz/2) 
        clear F1
        Ft  = fftshift(fft(f,[],1),1);     % fft(f) 
        F2t = fftshift(fft(f.*f,[],1),1);  % fft(f^2)   
        clear f
        % M(f0(z+backward_dz/2))  M1
        M1 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I+1)).*Ft,[],2),2),[],3),3);
        M1 = M1 + fftshift(fft(fftshift(fft(M_nonlinear(:,:,:,I+1).*F2t,[],2),2),[],3),3);
        clear   Ft  F2t
        F2 = excit_F.*expn2 + 0.25*backward_dz*expn2./(2i.*K).*(M + M1./expn2);% P12(z+backward_dz/2)
        clear M1
        F2(isnan(F2)) = 0;
        F2(real(Ktemp1)<=0) = 0;   
  
        % 3   
        f   = real(ifftn(ifftshift(F2))); % p1(z+backward_dz/2)   
        Ft  = fftshift(fft(f,[],1),1);    % fft(f1) 
        F2t = fftshift(fft(f.*f,[],1),1); % fft(f1^2)   
        clear f
        % M(f1(z+backward_dz/2)) M2
        M2 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I+1)).*Ft,[],2),2),[],3),3);
        M2 = M2 + fftshift(fft(fftshift(fft(M_nonlinear(:,:,:,I+1).*F2t,[],2),2),[],3),3);
        clear   Ft  F2t
        F3 = F2.*expn2 + 0.5*backward_dz*M2.*expn2./(2i.*K);  % P03(z+backward_dz)
        F3(isnan(F3)) = 0;
        F3(real(Ktemp1)<=0) = 0;  
    
        % 4     
        f   = real(ifftn(ifftshift(F3))); % p0(z+backward_dz)   
        Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
        F2t = fftshift(fft(f.*f,[],1),1); % fft(f2^2)
        clear f
        % M(f0(z+backward_dz))  M3
        M3 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I)).*Ft,[],2),2),[],3),3);
        M3 = M3 + fftshift(fft(fftshift(fft(M_nonlinear(:,:,:,I).*F2t,[],2),2),[],3),3);
        F4 = F2.*expn2 + 0.25*backward_dz.*(M2 + M3./expn2).*expn2./(2i.*K); % P14(z+backward_dz)
        clear M3   Ft  F2t
        F4(isnan(F4)) = 0;
        F4(real(Ktemp1)<=0) = 0;    
      
        % 5   
        f   = real(ifftn(ifftshift(F4))); % p0(z+backward_dz) 
        clear F4
        Ft  = fftshift(fft(f,[],1),1);    % fft(f2) 
        F2t = fftshift(fft(f.*f,[],1),1); % fft(f2^2)
        clear f
        % M(f1(z+backward_dz)) M4
        M4 = fftshift(fft(fftshift(fft((M_linear(:,:,:,I)).*Ft,[],2),2),[],3),3);
        M4 = M4 + fftshift(fft(fftshift(fft(M_nonlinear(:,:,:,I).*F2t,[],2),2),[],3),3);
        F5 = excit_F.*exp(1i*K*backward_dz) + ...
             backward_dz/6.0.*(M + 4*M2./expn2 + M4./exp(1i*K*backward_dz)).*...
             exp(1i*K*backward_dz)./(2i.*K); % P25(z+backward_dz)
       
        clear M   M2   M4  Ft  F2t
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0; 
         
%% add phase correction
        wxy2   = repmat(mgrid.w.', 1,mgrid.num_x, mgrid.num_y).*...
                 repmat(mgrid.w.', 1,mgrid.num_x, mgrid.num_y);
        e1     = exp(1i*(K - K2)*backward_dz);
        k_dif2 = (wxy2./(medium.c0*medium.c0) - wxy2./(cw(:,:,:,I).*cw(:,:,:,I)))*backward_dz;        
        k_dif  = (wxy2./(medium.c0*medium.c0) - wxy2./(cw(:,:,:,I+1).*cw(:,:,:,I+1)))*backward_dz;               
        
        correction = excit_F.*exp(1i*K2*backward_dz).*...
                    (1-e1.*(1+k_dif./4i./K)./(1-k_dif2./4i./K));
        clear k_dif  k_dif2  el  wxy2        
        
        F5 = F5 + correction;  
        clear correction
        F5(isnan(F5)) = 0;
        F5(real(Ktemp1)<=0) = 0;
        f = real(ifftn(ifftshift(F5)));
        
        rho1(1,:,:) = rho_rho(:,:,I+1);
        rho2(1,:,:) = rho_rho(:,:,I);   
        rho11 = repmat(rho1, mgrid.num_t, 1,1); 
        rho22 = repmat(rho2, mgrid.num_t, 1,1);            
        T_rhoc = 2.*rho11.*cw(:,:,:,I+1)./...
                (rho11.*cw(:,:,:,I+1) + rho22.*cw(:,:,:,I));          
        clear rho11 rho22 rho1          
        f = f./T_rhoc;
        
        F5 = fftshift(fftn(f));       
%%         
        excit_F = F5;
        f = real(ifftn(ifftshift(F5))); % p0(z+backward_dz) 
        clear F5
        
        if max(max(max(max(abs(f)))))==0
            disp ('The computation has lead to unphysical results. Please try reducing your step size.')
            break
        end 
        
        if reflection_order~=0
            Ktemp3 = (repmat(mgrid.w.',1,mgrid.num_x, mgrid.num_y)./...
                     cw(:,:,:,I+1)).^2 - ...
                     repmat(Ktemp1_kx.^2, mgrid.num_t, 1, mgrid.num_y) -...
                     repmat(Ktemp1_ky.^2, mgrid.num_t, mgrid.num_x, 1);                          
            K3 = sqrt(Ktemp3);
            if mod(mgrid.num_t,2) == 1
                K3(mgrid.num_t/2+0.5:end,:,:) = -K3(mgrid.num_t/2+0.5:end,:,:);
            else
            K3(mgrid.num_t/2:end,:,:)    = -K3(mgrid.num_t/2:end,:,:);
            end
            K3(K3==0) = eps;
            clear Ktemp3              
            
            P_in = excit_F.*exp(1i.*K3*mgrid.dz);
            p_ref(:,:,:, I) = real(ifftn(ifftshift(P_in))).*(T_rhoc-1);    
            clear P_in              
        end        
        
         if MDM_ani % show the animation
             p_ani(:,:,I) = squeeze(f(:,:,round(mgrid.num_y/2)));
             
             % in case a small computational time domain is used
             it = round(I*mgrid.dz/medium.c0/mgrid.dt);
             
             imagesc (mgrid.z*1e3, mgrid.x*1e3, squeeze(p_ani(it,:,:)))
             axis image
             drawnow;  
             if MDM_mov % record the animation
                 writeVideo(writerObj, getframe(gcf));
             end
         end  

        if sum(sum(sum(sensor_mask(:,:,I))))~=0
           % calculate the starting position of recorded signal
           sensor_mask_I1 = sum(sum(sum(sensor_mask(:,:,1:I-1)))) + 1;
    
           % calculate the ending position of recorded signal
           sensor_mask_I2 = sensor_mask_I1 -1 +  ...
           sum(sum(sum(sensor_mask(:,:,I))));
    
           % save the time-domain signals
           p_time(:,sensor_mask_I1:sensor_mask_I2) = f(:,sensor_mask(:,:,I)~=0);
        end
         
        clear p_I
         
        waitbar((mgrid.num_z - I+1)/mgrid.num_z)
     end
    close(h)    
    
end


%% flag for reflection
if (max(max(max(medium.c)))   ~= min(min(min(medium.c)))) || ...
    (max(max(max(medium.rho))) ~= min(min(min(medium.rho))))
    reflection_option = 1; 
else
    reflection_option = 0;
end

% calcuate the reflections
if reflection_order ~=0 && (reflection_option ==1)
    
    
    if MDMC == 1 % reflection with correction
        reflection = BMReflection3D(mgrid, medium, p_ref, ...
                     M_linear, K, cw, sensor_mask,...
                     reflection_order, excit_ps, M_nonlinear,...
                     MDMS, medium_cond);

    else  % reflection without correction
        reflection = BReflection3D(mgrid, medium, p_ref, ...
                     M_linear, K, cw, sensor_mask,...
                     reflection_order, excit_ps, M_nonlinear,...
                     MDMS, medium_cond);
    end
    
    if MDMS == 1  % with reflection correction
        p_time = reflection;  
    else    % without reflection correction
        p_time = p_time + reflection;  
    end        
           
end

end