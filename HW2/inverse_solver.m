
function[M]=HW2(error)

%x = [0 11 15 6 -7 3]';
%y = [0 0 6 13 10 -7]';
%d = [0.103 0.162  0.065  0.036 0.025 0.169]';

Xs = 8.1356;
Ys = -5.1398;
Zs = 11.5038;
P = 30.3326;
 
M=[Xs,Ys,Zs,P]';
%M=[3 -7 10 20]';

%lambda = 1e-5;
for ii = 1:1:1000
    
r=compute_residue(x,y,M,d);

%disp(norm(r));

[Grad,Hess]=compute_gradient_approx_hess(x,y,M,r);

%deltaM = (Hess+lambda*eye(4))\Grad';
deltaM= (Hess)\Grad';

M=M-deltaM;

if (norm(r)<1e-5)
    break;
end
%scatter(ii,norm(r));
end
end 



