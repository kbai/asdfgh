function [ Gamma, Hess] = compute_gradient_approx_hess( x,y,M,misfit )

xs = M(1);
ys = M(2);
zs = M(3);
p = M(4);

eta = ((x - xs).^2 + (y - ys).^2 + zs^2);

dx = x-xs;
dy = y-ys;

Ghat(:,1) = (3.*p.*zs.*(dx))./((eta).^(5/2));
Ghat(:,2) = (3.*p.*zs.*(dy))./((eta).^(5/2));
Ghat(:,3) = p./(eta).^(3/2) - (3*p.*zs.^2)./(eta).^(5/2);
Ghat(:,4) = zs./(eta).^(3/2);

sigma_P = 6;

P_mean = 35;

Gamma = (-(misfit')* Ghat)';

Gamma = Gamma + [0 0 0 1/p + (p-P_mean)/(sigma_P^2)]';

Hess = (Ghat')*Ghat ;

Hess(4,4) = Hess(4,4) - 1/p^2 + 1/(sigma_P^2) ;

%this is the apprximated Hessian;

end


