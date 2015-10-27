function [ Grad,Hess ,Hess_real] = compute_gradient_approx_hess( x,y,M,residue)

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


function [ r ] = compute_residue( x,y,M,d )
%this is calculating r=G(m)-d
xs = M(1);
ys = M(2);
zs = M(3);
p = M(4);
r = p*zs./((x - xs).^2 + (y - ys).^2 + zs^2).^(3/2)-d;
end


x = [0 11 15 6 -7 3]';
y = [0 0 6 13 10 -7]';
d = [0.103 0.162  0.065  0.036 0.025 0.169]';

Xs = 0;
Ys = 0;
Zs = 10;
P = 10;
M=[Xs,Ys,Zs,P]';

for ii = 1:1:100

r = compute_residue(x,y,M,d);

disp(norm(r));

[Grad,Hess]=compute_gradient_approx_hess(x,y,M,r);

deltaM= (Hess+1E-3*eye(4))\Grad';

M=M-deltaM;

end

