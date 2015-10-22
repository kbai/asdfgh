
function[M]=nonlinear_solver(x,y,d,Minit,w)

M=Minit;
r = 0;
r_old = 0;

for ii = 1:1:1000
r_old = r;    
r=compute_residue(x,y,M,d);


[Grad,Hess]=compute_gradient_approx_hess(x,y,M,r,w);

deltaM= (Hess)\Grad';

M=M+deltaM;

if (norm(r-r_old)<1e-7)
    break;
end
end
disp(r)
disp(ii)
end 



