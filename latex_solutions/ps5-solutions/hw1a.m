function hw2p1()
format long g
xi=[0 10 15 6 -7 3]';
yi=[0 0 6 13 10 7]';
zi=[0 0 0 0 0 0]';
ti=[322.418 321.031 321.228 323.093 324.415 322.706]';

% perturb the initial ti by a small amount
ti_new=[  322.42173921685
    321.033451805079
    321.231385791978
    323.082219348377
    324.407698376676
    322.696836727301];

% use a very good initial guess
M0=[315.147372203449  30.3124981621377 -17.1481612164449 15.9867180612486 5.26992727846599]';
[Mapp,errapp]=do_one_ti(xi,yi,zi,ti_new, M0,true, 100);
[M,err]      =do_one_ti(xi,yi,zi,ti_new, M0,false, 100);

% the exact one fail to converge to best result 
[Mapp, M]


% if we check the hessian after just the first step
[M,err,hess,hess_app,grad]      =do_one_ti(xi,yi,zi,ti_new, M0,false,1);
% exact hessian has a eigenvalue that's very small, so the exact hessian is
% near singular, the approximate hessian is better
[eigs(hess_app), eigs(hess)]
% a near singular hessian gives large update, which moves the solution away
% from the best guess, and converge to a local minimum
[hess_app\grad, hess\grad]


function [M,err,hess,hess_app,grad]=do_one_ti(xi,yi,zi,ti,M0,use_app,maxstep)

M=M0;

for step = 1:maxstep
    grad=zeros(5,1);
    hess=zeros(5,5);     % exact hessian
    hess_app=zeros(5,5); % approximate hessian
    for i=1:length(xi)        
        [tmp_grad, tmp_hess, tmp_hess_app]=get_grad_hessian(xi(i),yi(i),zi(i),ti(i),M);
        grad = grad + tmp_grad;
        hess = hess + tmp_hess;   
        hess_app = hess_app + tmp_hess_app;
    end        
    err(step) = norm(predict(xi,yi,zi,M)-ti);    
    if(step > 1 && err(step) > err(step-1))
        break;
    end
    if(use_app)
        M= M- pinv(hess_app)*grad;
    else
        M= M- pinv(hess)*grad;
    end
end

function [grad,hess,hess_approx] = get_grad_hessian(xi,yi,zi,ti,M)
ts=M(1); xs=M(2); ys=M(3); zs=M(4); v=M(5);
R=sqrt((xs-xi)^2 + (ys-yi)^2 + (zs-zi)^2);
nx=(xs-xi)/R;
ny=(ys-yi)/R;
nz=(zs-zi)/R;
e=(ts+R/v-ti);

grad = e*[1 nx/v ny/v nz/v -R/v^2]';
sen = [1 nx/v ny/v nz/v -R/v^2]';
hess_approx = sen*sen'; % approximate hess

hess=[1 nx/v ny/v nz/v -R/v^2
    0  1/v*(1/v*nx^2+e/R*(1-nx^2))  1/v*(1/v*nx*ny-e/R*nx*ny) 1/v*(1/v*nx*nz-e/R*nx*nz) -1*nx/v^2*(R/v+e)
    0  0  1/v*(1/v*ny^2+e/R*(1-ny^2))  1/v*(1/v*ny*nz-e/R*ny*nz) -1*ny/v^2*(R/v+e)
    0  0   0  1/v*(1/v*nz^2+e/R*(1-nz^2))  -1*nz/v^2*(R/v+e)
    0 0 0 0 1*R/v^3*(R/v+2*e)];
hess = hess + hess';
for k=1:5
    hess(k,k)=hess(k,k)/2;
end

function r=predict(xi,yi,zi,M)
ts=M(1); xs=M(2); ys=M(3); zs=M(4); v=M(5);
r=xi*0;
for i=1:length(xi)
    R=sqrt((xs-xi(i))^2 + (ys-yi(i))^2 + (zs-zi(i))^2);
    r(i)=(ts+R/v);
end
    