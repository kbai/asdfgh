function hw2_1
set(0,'defaulttextfontname','times','defaulttextfontsize',14);
set(0,'defaultaxesfontname','times','defaultaxesfontsize',14);

sigma = 0.001;
% std_mc = [0.099670 0.137469 0.149719 0.691071]';
%%%hw 2 problem 1d
xi = [0 11 15 6 -7 3]';
yi = [0 0 6 13 10 -7]';
ui = [0.103 0.162  0.065  0.036 0.025 0.169]';
M0 = [8 -5 10 30]';  %initial guess, xs ys zs P

[Ms, mcov] = nonlinear_solver(xi,yi,ui,M0,sigma); 

% model covariance matrix
disp('model covariance matrix:');
disp(mcov);
disp('standard deviation');
std_m = sqrt(diag(mcov))

% calculate correlation matrix
s = diag(std_m);
disp('coefficient matrix:');
mcoef = s^(-1)*mcov*s^(-1)

end

function[M, mcov]=nonlinear_solver(x,y,ui,Minit,sigma)

M = Minit;

misfit_old = 0;

misfit = 0;

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

mcov = inv(Hess/sigma^2);

%disp(['Number of iterations:',num2str(ii)]);
end

function [ misfit ] = compute_misfit( x,y,M,ui )

xs = M(1);
ys = M(2);
zs = M(3);
p = M(4);

misfit  = ui - p*zs./((x - xs).^2 + (y - ys).^2 + zs^2).^(3/2);

end

function [ Grad, Hess] = compute_gradient_approx_hess( x,y,M,misfit )

xs = M(1);
ys = M(2);
zs = M(3);
p = M(4);

eta = ((x - xs).^2 + (y - ys).^2 + zs^2);

dx = x-xs;
dy = y-ys;

Ghat(:,1) = (3.*p.*zs.*(dx))./((eta).^(5/2));
Ghat(:,2) = (3.*p.*zs.*(dy))./((eta).^(5/2));
Ghat(:,3) = p./(eta).^(3/2) - (3*p.*zs.^2)./(eta).^(5/2);
Ghat(:,4) = zs./(eta).^(3/2);

Grad = (misfit')* Ghat;

Hess = (Ghat')*Ghat;
%this is the apprximated Hessian;

end
