function [ Grad,Hess ] = compute_gradient_approx_hess( x,y,M,residue)
%COMPUTE_GRADIENT_APPROX_HESS Summary of this function goes here
%   Detailed explanation goes here
xs = M(1);
ys = M(2);
zs = M(3);
p = M(4);

R = ((x - xs).^2 + (y - ys).^2 + zs^2);

Jacob(:,1) = (3.*p.*zs.*(2*x - 2*xs))./(2*(R).^(5/2));
Jacob(:,2) = (3.*p.*zs.*(2*y - 2*ys))./(2*(R).^(5/2));
Jacob(:,3) = p./(R).^(3/2) - (3*p.*zs.^2)./(R).^(5/2);
Jacob(:,4) = zs./(R).^(3/2);
disp('finished');

Grad = 2*(residue')* Jacob;
%this is the apprximated Hessian;
Hess = 2*(Jacob')*Jacob;
Hess = 0.5*(Hess + Hess');

end

