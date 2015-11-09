function[con1,con2,M]=nonlinear_solver(x,y,ui,Minit)

M = Minit;
misfit_old = 0;
misfit_r = zeros(1,100);
Ms_r1 = zeros(4,100);
Ms_r2 = zeros(4,100);
n_iter1 = 0;
n_iter2 = 0;

for ii = 1:1:1000
    
    misfit = compute_misfit(x,y,M,ui);
    
    [Gamma,Hess] = compute_gradient_approx_hess2(x,y,M,misfit);
    
    deltaM = - (Hess)\Gamma;
    
    Ms_r1(:,ii) = M;
    
    M=M + deltaM;
    
    if (norm(misfit - misfit_old)<1e-7)
        
        Ms_r1(:,ii+1) = M;

        break;
    end
    misfit_r(ii) = misfit'*misfit;
    

    
    misfit_old = misfit;
end

disp(['Number of iterations:',num2str(ii)]);
n_iter1 = ii;
plot(1:1:100,log(misfit_r)/log(10));
hold on
con1=(norm(misfit)<0.02);
M = Minit;
misfit_old = 0;
misfit_r = zeros(1,100);
for ii = 1:1:1000
    
    misfit = compute_misfit(x,y,M,ui);
    
    [Gamma,Hess] = compute_gradient_approx_hess(x,y,M,misfit);
    
    deltaM = - (Hess)\Gamma;
    
    Ms_r2(:,ii) = M;

    M=M + deltaM;
    
    if (norm(misfit - misfit_old)<1e-7)
        Ms_r2(:,ii+1) = M;

        
        break;
    end
    misfit_r(ii) = misfit'*misfit;
    
    misfit_old = misfit;

end

n_iter2 = ii;
disp(['Number of iterations:',num2str(ii)]);

plot(1:1:100,log(misfit_r)/log(10));
legend('exact hessian','approximate hessian')
xlabel('number of iterations');
ylabel('log_{10}(misfit^2)');
hold off

figure(2)
plot(Ms_r1(1,1:n_iter1),Ms_r1(4,1:n_iter1),'o-');

hold on
plot(Ms_r2(1,1:n_iter2),Ms_r2(4,1:n_iter2),'o-');

con2=(norm(misfit)<0.002);
con2 = 2*(n_iter2>n_iter1)+(n_iter1==n_iter2);

end




