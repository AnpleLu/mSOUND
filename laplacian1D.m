function v = laplacian1D(f, dz)

% DESCRIPTION:
%  the centered difference formulas for five-point stencils approximating
%  the second order derivatives

% USAGE:
% v = laplacian1D(f, dz)

% INPUTS:
% f        input 1D matrix f
% dz       step size

% OUTPUTS:
% v       second order derivative of f

%%
[n m] = size(f);
vz = zeros(size(f));

%2nd order on the boundary
vz(1,:) = (2*f(1,:)-5*f(2,:)+4*f(3,:)-f(4,:))/dz^2;
vz(2,:) = (f(1,:)-2*f(2,:)+f(3,:))/dz^2;
vz(n-1,:) = (f(n-2,:)-2*f(n-1,:)+f(n,:))/dz^2;
vz(n,:) = (2*f(n,:)-5*f(n-1,:)+4*f(n-2,:)-f(n-3,:))/dz^2;

for i=3:n-2
    
    vz(i,:) = (-f(i+2,:)+16*f(i+1,:)-30*f(i,:)+16*f(i-1,:)-f(i-2,:))/12/dz^2;

end


v = vz;


