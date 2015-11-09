function [ Gamma, Hess] = compute_gradient_approx_hess2( x,y,M,misfit)

xs = M(1);
ys = M(2);
zs = M(3);
p = M(4);

R = ((x - xs).^2 + (y - ys).^2 + zs^2);

dx = x-xs;
dy = y-ys;

Ghat(:,1) = (3.*p.*zs.*(dx))./((R).^(5/2));
Ghat(:,2) = (3.*p.*zs.*(dy))./((R).^(5/2));
Ghat(:,3) = p./(R).^(3/2) - (3*p.*zs.^2)./(R).^(5/2);
Ghat(:,4) = zs./(R).^(3/2);

Gamma = - Ghat'*misfit;
%this is the apprximated Hessian;


dHess = zeros(4);
dHess(1,2) = - misfit'*((15*p*zs*(dx).*(dy))./((R).^(7/2)));
dHess(1,3) = - misfit'*((3*p.*(dx))./((R).^(5/2)) - (15*p*zs.^2.*(dx))./((R).^(7/2)));
dHess(1,4) = - misfit'*((3*zs.*(dx))./((R).^(5/2)));
dHess(2,3) = - misfit'*( (3*p.*(dy))./((R).^(5/2)) - (15*p*zs.^2.*(dy))./((R).^(7/2)));
dHess(2,4) = - misfit'*((3*zs.*(dy))./((R).^(5/2)));
dHess(3,4) = - misfit'*(1./(R).^(3/2) - (3*zs^2)./(R).^(5/2));

dHess = (dHess + dHess');

dHess(1,1) = - misfit'*((15*p*zs.*(dx).^2)./((R).^(7/2)) - (3*p.*zs)./(R).^(5/2));
dHess(2,2) = - misfit'*((15*p*zs.*(dy).^2)./((R).^(7/2)) - (3*p.*zs)./(R).^(5/2));
dHess(3,3) = - misfit'*((15*p*zs^3)./(R).^(7/2) - (9*p*zs)./(R).^(5/2));
dHess(4,4) = 0;



Hess = (Ghat')*Ghat;
%only add dHess if using exact Hessian
Hess = Hess + dHess;


end


