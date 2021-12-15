function [medium_cond] = medium_case(medium)
% DESCRIPTION:
% decide the medium properties

% USAGE:
% [medium_case] = medium_case(medium)

% INPUTS:
% mgrid        structure to define the computation domain
% medium       medium properties

% OUTPUTS:
% medium_cond 
%%
if max(max(max(medium.c)))   == min(min(min(medium.c)))&& ... %speed and density are constants
   max(max(max(medium.rho))) == min(min(min(medium.rho)))&&...
   max(max(max(medium.ca)))  == min(min(min(medium.ca)))&&... 
   max(max(max(medium.cb)))  == 2.0&&... 
   min(min(min(medium.cb)))  == 2.0&&...
   max(max(max(medium.beta)))== 0 % Mnonlnear = 0  
medium_cond = 1; % 'linear homogeneous medium'

elseif max(max(max(medium.c)))   == min(min(min(medium.c)))&& ...% Mlinear= c
       max(max(max(medium.rho))) == min(min(min(medium.rho)))&&... 
       max(max(max(medium.ca)))  == min(min(min(medium.ca)))&&...
       max(max(max(medium.cb)))  == 2.0&&...
       min(min(min(medium.cb)))  == 2.0&&...
       max(max(max(medium.beta)))~= 0&&... % Mnonlnear = c
       max(max(max(medium.beta)))== min(min(min(medium.beta)))    
medium_cond = 2; % 'nonlinear homogeneous medium'

elseif max(max(max(medium.c)))   == min(min(min(medium.c)))&& ...% Mlinear = 0
       max(max(max(medium.rho))) == min(min(min(medium.rho)))&&...
       max(max(max(medium.ca)))  == min(min(min(medium.ca)))&&... 
       max(max(max(medium.cb)))  == min(min(min(medium.cb)))&&....
       min(min(min(medium.cb)))  == 2.0&&...
       max(max(max(medium.beta)))~= min(min(min(medium.beta))) % Mnonlnear ~= 0
medium_cond = 3; % 'only nonlinear coefficient variation'

elseif (max(max(max(medium.c)))   ~= min(min(min(medium.c)))|| ...% Mlinear~=0
        max(max(max(medium.rho))) ~= min(min(min(medium.rho)))||...
        max(max(max(medium.ca)))  ~= min(min(min(medium.ca))))&&...
        max(max(max(medium.beta)))== 0 % Mnonlnear = 0   
medium_cond = 4; % 'linear media'
   
else    
medium_cond = 5; % 'other media'

end




end