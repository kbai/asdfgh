function [ Grad, Hess] = compute_gradient_approx_hess2( x,y,M,residue)

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
%this is the apprximated Hessian;


%dHess = (G(m)-d)^T*Q
dHess = zeros(4);
dHess(1,2) = residue'*((15*p*zs*(dx).*(dy))./(4*(R).^(7/2)));
dHess(1,3) = residue'*((3*p.*(dx))./(2*(R).^(5/2)) - (15*p*zs.^2.*(dx))./(2*(R).^(7/2)));
dHess(1,4) = residue'*((3*zs.*(dx))./(2*(R).^(5/2)));
dHess(2,3) = residue'*( (3*p.*(dy))./(2*(R).^(5/2)) - (15*p*zs.^2.*(dy))./(2*(R).^(7/2)));
dHess(2,4) = residue'*((3*zs.*(dy))./(2*(R).^(5/2)));
dHess(3,4) = residue'*(1./(R).^(3/2) - (3*zs^2)./(R).^(5/2));

dHess = (dHess + dHess');

dHess(1,1) = residue'*((15*p*zs.*(dx).^2)./(4*(R).^(7/2)) - (3*p.*zs)./(R).^(5/2));
dHess(2,2) = residue'*((15*p*zs.*(dy).^2)./(4*(R).^(7/2)) - (3*p.*zs)./(R).^(5/2));
dHess(3,3) = residue'*((15*p*zs^3)./(R).^(7/2) - (9*p*zs)./(R).^(5/2));
dHess(4,4) = 0;



Hess = (Jacob')*Jacob;
%only add dHess if using exact Hessian
Hess = Hess+dHess;

Hess = 0.5*(Hess + Hess');

end


