syms x y p xs ys zs d;
u = (p*zs/((x-xs)^2+(y-ys)^2+zs^2).^(3/2)-d)^2;
f1=diff(u,xs);
f2=diff(u,ys);
f3=diff(u,zs);
f4=diff(u,p);

g11=diff(diff(u,xs),xs);
g12=diff(diff(u,xs),ys);

g13=diff(diff(u,xs),zs);

g22=diff(diff(u,ys),ys);
g23=diff(diff(u,ys),zs);

g33=diff(diff(u,zs),zs);






latex(f1)
latex(f2)
latex(f3)
latex(f4)

x = [0 11 15 6 -7 3]';
y = [0 0 6 13 10 -7]';
d = [0.103 0.162  0.065  0.036 0.025 0.169]';

Xs = 0;
Ys = 0;
Zs = 10;
P = 10;
M=[Xs,Ys,Zs,P]';

for ii = 1:1:10
r=compute_residue(x,y,M,d);
disp(norm(r));
[Grad,Hess]=compute_gradient_approx_hess(x,y,M,r);
deltaM= (Hess+1E-3*eye(4))\Grad';
M=M-deltaM;
end

