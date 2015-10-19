
function[M]=nonlinear_solver(x,y,d,Minit)

%x = [0 11 15 6 -7 3]';
%y = [0 0 6 13 10 -7]';
%d = [0.103 0.162  0.065  0.036 0.025 0.169]'+error;
M=Minit;
%M=[3 -7 10 20]';

%lambda = 1e-5;
for ii = 1:1:1000
    
r=compute_residue(x,y,M,d);

%disp(norm(r));

[Grad,Hess]=compute_gradient_approx_hess(x,y,M,r);

%deltaM = (Hess+lambda*eye(4))\Grad';
deltaM= (Hess)\Grad';

M=M-deltaM;

if (norm(r)<1e-7)
    break;
end
end
end 



