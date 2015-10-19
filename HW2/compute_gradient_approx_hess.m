function [ Grad, Hess] = compute_gradient_approx_hess( x,y,M,residue)

xs = M(1);
ys = M(2);
zs = M(3);
p = M(4);

R = ((x - xs).^2 + (y - ys).^2 + zs^2);

dx = x-xs;
dy = y-ys;

Jacob(:,1) = (3.*p.*zs.*(dx))./((R).^(5/2));
Jacob(:,2) = (3.*p.*zs.*(dy))./((R).^(5/2));
Jacob(:,3) = p./(R).^(3/2) - (3*p.*zs.^2)./(R).^(5/2);
Jacob(:,4) = zs./(R).^(3/2);

Grad = (residue')* Jacob;



Hess = (Jacob')*Jacob;
%this is the apprximated Hessian;


Hess = 0.5*(Hess + Hess');

end


