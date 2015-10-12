syms x y p xs ys zs;
u = p*zs/((x-xs)^2+(y-ys)^2+zs^2).^(3/2);
diff(u,xs);
diff(u,ys);
diff(u,zs);
diff(u,p);

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
deltaM= Hess\Grad';
M=M-deltaM;
end

