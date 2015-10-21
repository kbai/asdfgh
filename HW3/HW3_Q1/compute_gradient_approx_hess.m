function [ Grad, Hess] = compute_gradient_approx_hess( x,y,M,residue,weigh)

W =diag(weigh);

C=W'*W;

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

Grad = (residue')*C*Ghat;



Hess = (Ghat')*C*Ghat;
%this is the apprximated Hessian;


Hess = 0.5*(Hess + Hess');

end


