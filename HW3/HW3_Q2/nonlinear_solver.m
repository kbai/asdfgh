
function[M]=nonlinear_solver(x,y,ui,Minit,w)

M=Minit;
misfit = 0;
misfit_old = 0;

for ii = 1:1:1000
misfit_old = misfit;    
misfit=compute_misfit(x,y,M,ui);


[Grad,Hess]=compute_gradient_approx_hess(x,y,M,misfit,w);

deltaM= (Hess)\Grad';

M=M+deltaM;

if ((misfit-misfit_old)'*(misfit-misfit_old)<1e-7)
    break;
end
end
disp(misfit)
disp(ii)
end 



