% this starts from 03/03/2017;
% it plots the structure with three different cylinders;
%% nodes per wavelength
function[c, rho, beta, ca] = annulus(dx, dz, Nx, Nz, ...
                                fc, c0,...
                                speed_water, speed_fat,...
                                rho_water, rho_fat,...
                                beta_water, beta_fat,...
                                ca_water, ca_fat, p_shift)
                 
lambda = c0/fc;
[x_mesh,z_mesh] = ndgrid((([1:Nx]-Nx/2-1/2)*dx), ...
			             (([1:Nz]-Nz/2)*dz));

%%   largest skull model
% r1c = 41*lambda;
% r2c = 28*lambda;                    
% r1 = sqrt((x_mesh-0*lambda).^2 + (z_mesh - 32*lambda-p_shift).^2);
% r2 = sqrt((x_mesh-0*lambda).^2 + (z_mesh - 32*lambda-p_shift).^2);

%%   middle size skull model
r1c = 38*lambda;
r2c = 32*lambda;                    
r1 = sqrt((x_mesh-0*lambda).^2 + (z_mesh - 32*lambda-p_shift).^2);
r2 = sqrt((x_mesh-0*lambda).^2 + (z_mesh - 32*lambda-p_shift).^2);

%%   small size skull model
% r1c = 35*lambda;
% r2c = 32*lambda;                    
% r1 = sqrt((x_mesh-0*lambda).^2 + (z_mesh - 32*lambda-p_shift).^2);
% r2 = sqrt((x_mesh-0*lambda).^2 + (z_mesh - 32*lambda-p_shift).^2);

%%   smallest size skull model
% r1c = 37*lambda;
% r2c = 35.5*lambda;                    
% r1 = sqrt((x_mesh-0*lambda).^2 + (z_mesh - 32*lambda-p_shift).^2);
% r2 = sqrt((x_mesh-0*lambda).^2 + (z_mesh - 32*lambda-p_shift).^2);

%%
c = speed_water*(r1>r1c | r2<=r2c)+...
    speed_fat*(r1<=r1c & r2>r2c);

rho = rho_water*(r1>r1c | r2<=r2c)+...
      rho_fat*(r1<=r1c & r2>r2c);
  
beta = beta_water*(r1>r1c | r2<=r2c)+...
       beta_fat*(r1<=r1c & r2>r2c);
 
ca = ca_water*(r1>r1c | r2<=r2c)+...
     ca_fat*(r1<=r1c & r2>r2c);
 
% gradient along the boundary of the irregularty

end



