function [excit_p] = excitation(excit_ps, mgrid)
% DESCRIPTION:
%     format the user-defined excitation for simulation implementation

% USAGE:
%     [excit_p] = excitation(excit_ps, mgrid)
%
% INPUTS:
% excit_ps  user-defined time domain signal
% mgrid

% OUTPUTS:
% excit_p   formatted excitation signal in the time domain
%%
% 1D simulation
if mgrid.dim == 1 
    
    % excitation pulse length
    [length_ts] = length(excit_ps);
    if length_ts ~= 1 % transient simulation
        excit_p = [excit_ps; zeros(mgrid.num_t-length_ts, 1)];
    else              % frequency-domain simulation
        excit_p = excit_ps;
    end
    
% 2D simulation    
elseif mgrid.dim == 2
    
    % excitational spulse size
    [length_ts,num_elx] = size(excit_ps);
    if length_ts ~= 1 % transient simulation
        excit_p = [excit_ps; zeros(mgrid.num_t-length_ts, num_elx)];  
        excit_ptemp = zeros(mgrid.num_t, mgrid.num_x);
        excit_ptemp(:,round((mgrid.num_x-num_elx)/2)+1:...
                      round((mgrid.num_x-num_elx)/2)+num_elx) = excit_p;  
    else  % frequency-domain simulation
        excit_ptemp = zeros(	1, mgrid.num_x);
        excit_ptemp(:,round((mgrid.num_x-num_elx)/2)+1:...
                      round((mgrid.num_x-num_elx)/2)+num_elx) = excit_ps;     
    end
    excit_p = excit_ptemp;

% 3D simulation
elseif mgrid.dim == 3
    
    % excitation size
    [length_ts, num_elx, num_ely] = size(excit_ps);
    
    if length_ts~=1 % transient simulation
        
        excit_p = [excit_ps; zeros(mgrid.num_t-length_ts, num_elx, num_ely)];
        excit_ptemp = zeros(mgrid.num_t, mgrid.num_x, mgrid.num_y);
        excit_ptemp(:,round((mgrid.num_x-num_elx)/2)+1:...
                      round((mgrid.num_x-num_elx)/2)+num_elx,...
                      round((mgrid.num_y-num_ely)/2)+1:...
                      round((mgrid.num_y-num_ely)/2)+num_ely) = excit_p;     
    else   % frequency-domain simulation
        excit_ptemp = zeros(1, mgrid.num_x, mgrid.num_y);
        excit_ptemp(:,round((mgrid.num_x-num_elx)/2)+1:...
                      round((mgrid.num_x-num_elx)/2)+num_elx,...
                      round((mgrid.num_y-num_ely)/2)+1:...
                      round((mgrid.num_y-num_ely)/2)+num_ely) = excit_ps;  
    
    end
    excit_p = excit_ptemp;

end

end
