
function[M]=nonlinear_solver(x,y,ui,Minit)

%% for damped least square method

Mr1 = zeros(1000,length(Minit));

M = Minit;

misfit_old = 0;

misfit = 0;

for ii = 1:1:20
    
misfit = compute_misfit(x,y,M,ui);


[Gamma,Hess] = compute_gradient_approx_hess(x,y,M,misfit);

deltaM= -(Hess+(ii<=3)*10*diag(diag(Hess)))\Gamma; 

Mr1(ii,:)=M;

M = M + deltaM;


if (norm(misfit - misfit_old)<1e-7)
    break;
end

misfit_old = misfit;

end

n1=ii;
%% for undamped least square method

Mr2 = zeros(1000,length(Minit));

M = Minit;

misfit_old = 0;

misfit = 0;

for ii = 1:1:2
 %% since LS method do not converge, we plot only the first 2 points.
 
misfit = compute_misfit(x,y,M,ui);


[Gamma,Hess] = compute_gradient_approx_hess(x,y,M,misfit);

deltaM= -(Hess)\Gamma; 

Mr2(ii,:)=M;

M = M + deltaM;


if (norm(misfit - misfit_old)<1e-7)
    break;
end

misfit_old = misfit;

end
n2=ii;

disp(['Number of iterations:',num2str(ii)]);
%% plot the contour and convergence path
[Xs,P]=meshgrid(7.5:0.01:10.5,5:0.1:35);
[m,n]=size(Xs);
ERROR = Xs*0;
for ll = 1:1:m
    for kk = 1:1:n
        M0=[Xs(ll,kk),-5.1412,11.5066,P(ll,kk)]';
   %     M0 = [Xs(ll,kk),-10,10,P(ll,kk)];
        Misfit = compute_misfit(x,y,M0,ui);
        ERROR(ll,kk) = Misfit'*Misfit;
    end
end
contour(Xs,P,ERROR,300);
hold on
plot(Mr1(1:n1,1),Mr1(1:n1,4),'-o','linewidth',2);
hold on
plot(Mr2(1:n2,1),Mr2(1:n2,4),'-o','linewidth',2);
colorbar
legend('Contour','LM method','GN method')
print('fig1b.pdf','-dpdf');
end 



