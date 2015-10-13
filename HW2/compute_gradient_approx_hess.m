function [ Grad,Hess ,Hess_real] = compute_gradient_approx_hess( x,y,M,residue)
%COMPUTE_GRADIENT_APPROX_HESS Summary of this function goes here
%   Detailed explanation goes here
xs = M(1);
ys = M(2);
zs = M(3);
p = M(4);

R = ((x - xs).^2 + (y - ys).^2 + zs^2);

dx = x-xs;
dy = y-ys;

Jacob(:,1) = (3.*p.*zs.*(dx))./(2*(R).^(5/2));
Jacob(:,2) = (3.*p.*zs.*(dy))./(2*(R).^(5/2));
Jacob(:,3) = p./(R).^(3/2) - (3*p.*zs.^2)./(R).^(5/2);
Jacob(:,4) = zs./(R).^(3/2);
disp('finished');

Grad = 2*(residue')* Jacob;
%this is the apprximated Hessian;



Hess = zeros(4);
Hess(1,2) =2 * residue'*((15*p*zs*(dx).*(dy))./(4*(R).^(7/2)));
Hess(1,3) =2 * residue'*((3*p.*(dx))./(2*(R).^(5/2)) - (15*p*zs.^2.*(dx))./(2*(R).^(7/2)));
Hess(1,4) =2 * residue'*((3*zs.*(dx))./(2*(R).^(5/2)));
Hess(2,3) =2 * residue'*( (3*p.*(dy))./(2*(R).^(5/2)) - (15*p*zs.^2.*(dy))./(2*(R).^(7/2)));
Hess(2,4) =2 * residue'*((3*zs.*(dy))./(2*(R).^(5/2)));
Hess(3,4) =2 * residue'*(1./(R).^(3/2) - (3*zs^2)./(R).^(5/2));

Hess = (Hess + Hess');

Hess(1,1) =2 * residue'*((15*p*zs.*(dx).^2)./(4*(R).^(7/2)) - (3*p.*zs)./(R).^(5/2));
Hess(2,2) =2 * residue'*((15*p*zs.*(dy).^2)./(4*(R).^(7/2)) - (3*p.*zs)./(R).^(5/2));
Hess(3,3) =2 * residue'*((15*p*zs^3)./(R).^(7/2) - (9*p*zs)./(R).^(5/2));
Hess(4,4) = 0;


Hess_real = Hess;


Hess = 2*(Jacob')*Jacob;
Hess = 0.5*(Hess + Hess');

end


