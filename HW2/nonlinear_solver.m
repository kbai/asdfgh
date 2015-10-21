
function[M]=nonlinear_solver(x,y,d,Minit)

M=Minit;

for ii = 1:1:1000
    
r=compute_residue(x,y,M,d);


[Grad,Hess]=compute_gradient_approx_hess(x,y,M,r);

deltaM= (Hess)\Grad';

M=M+deltaM;

if (norm(r)<1e-7)
    break;
end
end
end 



