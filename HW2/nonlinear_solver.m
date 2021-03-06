function[M, mcov]=nonlinear_solver(x,y,ui,Minit)

M = Minit;
misfit_old = 0;

for ii = 1:1:1000
    
    misfit = compute_misfit(x,y,M,ui);
    
    [Grad,Hess] = compute_gradient_approx_hess(x,y,M,misfit);
    
    deltaM= (Hess)\Grad';
    
    M=M + deltaM;
    
    if (norm(misfit - misfit_old)<1e-7)
        break;
    end
    
    misfit_old = misfit;
end

%mcov = inv(Hess/sigma^2);

%disp(['Number of iterations:',num2str(ii)]);
end



