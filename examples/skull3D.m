% this starts from 03/03/2017;
% it plots the structure with three different cylinders;
%% nodes per wavelength
function[c, rho, beta, ca] = skull3D(dx, dy, dz, Nx, Ny, Nz, ...
                                fc, c0,...
                                speed_water, speed_fat,...
                                rho_water, rho_fat,...
                                beta_water, beta_fat,...
                                ca_water, ca_fat, p_shift)
                 
lambda = c0/fc;
[x_mesh, y_mesh, z_mesh] = ndgrid((([1:Nx]-Nx/2-1/2)*dx), ...
			                     (([1:Ny]-Ny/2-1/2)*dy), ...
                                 (([1:Nz]-Nz/2)*dz));

                            
r1 = sqrt((x_mesh -4*lambda).^2 + (y_mesh -4*lambda).^2 + (z_mesh - 20*lambda-p_shift).^2);
r2 = sqrt((x_mesh- 4*lambda).^2 + (y_mesh -4*lambda).^2 + (z_mesh - 20*lambda-p_shift).^2);

radius1 = 23*lambda;
radius2 = 21*lambda;
c = speed_water*(r1>radius1 | r2<=radius2)+...
    speed_fat*(r1<=radius1 & r2>radius2);


rho = rho_water*(r1>radius1 | r2<=radius2)+...
      rho_fat*(r1<=radius1 & r2>radius2);
  
beta = beta_water*(r1>radius1 | r2<=radius2)+...
       beta_fat*(r1<=radius1 & r2>radius2);
 
ca = ca_water*(r1>radius1 | r2<=radius2)+...
     ca_fat*(r1<=radius1 & r2>radius2);
 
% gradient along the boundary of the irregularty

end



