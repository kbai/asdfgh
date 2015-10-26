
function[M]=nonlinear_solver(x,y,d,Minit)

M = Minit;

misfit_old = 0;

misfit = 0;

for ii = 1:1:1000
    
misfit = compute_misfit(x,y,M,d);


[Grad,Hess] = compute_gradient_approx_hess(x,y,M,misfit);

deltaM= (Hess)\Grad'; 

M=M + deltaM;

if (norm(misfit - misfit_old)<1e-7)
    break;
end

misfit_old = misfit;
end

disp(['Number of iterations:',num2str(ii)]);
end 



