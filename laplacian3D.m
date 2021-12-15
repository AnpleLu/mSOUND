function v = laplacian3D(f, dx, dy, dz)

% DESCRIPTION:
%  the centered difference formulas for five-point stencils approximating
%  the second order derivatives

% USAGE:
% v=laplacian3D(f, dx, dy, dz)

% INPUTS:
% f        input 3D matrix f
% dx       step size
% dy       step size
% dz       step size

% OUTPUTS:
% v       second order derivative of f

[n m p] = size(f);
vx = zeros(size(f));
vy = zeros(size(f));
vz = zeros(size(f));

%2nd derivative with regard to x
%2nd order on the boundary
vx(1,:,:) = (2*f(1,:,:)-4*f(2,:,:)+4*f(3,:,:)-f(4,:,:)-f(2,:,:))/dx^2;
vx(2,:,:) = (f(1,:,:)-2*f(2,:,:)+f(3,:,:))/dx^2;
vx(n-1,:,:) = (f(n-2,:,:)-2*f(n-1,:,:)+f(n,:,:))/dx^2;
vx(n,:,:)   = (2*f(n,:,:)-5*f(n-1,:,:)+4*f(n-2,:,:)-f(n-3,:,:))/dx^2;

for i=3:n-2
    vx(i,:,:) = (-f(i+2,:,:)+16*f(i+1,:,:)-30*f(i,:,:)+16*f(i-1,:,:)-f(i-2,:,:))/12/dx^2;
end

%2nd derivative with regard to y
%2nd order on the boundary
vy(:,1,:) = (2*f(:,1,:)-4*f(:,2,:)+4*f(:,3,:)-f(:,4,:) - f(:,2,:))/dy^2; % write in this way to reduce round-off error
vy(:,2,:) = (f(:,1,:)-2*f(:,2,:)+f(:,3,:))/dy^2;
vy(:,m-1,:) = (f(:,m-2,:)-2*f(:,m-1,:)+f(:,m,:))/dy^2;
vy(:,m,:) = (2*f(:,m,:)-5*f(:,m-1,:)+4*f(:,m-2,:)-f(:,m-3,:))/dy^2;
for i=3:m-2     
    vy(:,i,:) = (-f(:,i+2,:)+16*f(:,i+1,:)-30*f(:,i,:)+16*f(:,i-1,:)-f(:,i-2,:))/12/dy^2;
end

%2nd derivative with regard to z
%2nd order on the boundary
vz(:,:,1) = (2*f(:,:,1)-4*f(:,:,2)+4*f(:,:,3)-f(:,:,4) - f(:,:,2))/dz^2; % write in this way to reduce round-off error
vz(:,:,2) = (f(:,:,1)-2*f(:,:,2)+f(:,:,3))/dz^2;
vz(:,:,p-1) = (f(:,:,p-2)-2*f(:,:,p-1)+f(:,:,p))/dz^2;
vz(:,:,p) = (2*f(:,:,p)-5*f(:,:,p-1)+4*f(:,:,p-2)-f(:,:,p-3))/dz^2;

for i=3:p-2 
    vz(:,:,i) = (-f(:,:,i+2)+16*f(:,:,i+1)-30*f(:,:,i)+16*f(:,:,i-1)-f(:,:,i-2))/12/dz^2;
end


v = vx+vy+vz;

