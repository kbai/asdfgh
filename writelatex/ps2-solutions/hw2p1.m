function hw2p1()
rng(0);
ti=[322.418 321.031 321.228 323.093 324.415 322.706]';

for it=1:1000
    disp(it);
    ti_new = ti + 0.01*randn(6,1);
    
    % inexact Hessian
    [Mapp(:,it),tmp_err1]=do_one_ti(ti_new, true);
    err_app(it)=min(tmp_err1); 
    
    % exact Hessian
    [M(:,it),tmp_err2]=do_one_ti(ti_new, false);
    err(it)=min(tmp_err2);
end
scatter(Mapp(2,:), Mapp(3,:), 10, Mapp(1,:)); colorbar;
print -depsc hwp1_fig1_app.eps
hold on;
plot([0 10 15 6 -7 3],[0 0 6 13 10 7],'kv');
hold off;
print -depsc hwp1_fig2_app.eps

scatter(M(2,:), M(3,:), 10, M(1,:)); colorbar;
print -depsc hwp1_fig1.eps
hold on;
plot([0 10 15 6 -7 3],[0 0 6 13 10 7],'kv');
hold off;
print -depsc hwp1_fig2.eps

function [M,err]=do_one_ti(ti, use_app)
xi=[0 10 15 6 -7 3]';
yi=[0 0 6 13 10 7]';
zi=[0 0 0 0 0 0]';

%M=[315 30 -17 16 5]';

M=[315.147372203449           30.3124981621377          -17.1481612164449           15.9867180612486           5.26992727846599]';

M2=M;

for step = 1:1000
    grad=zeros(5,1);
    hess=zeros(5,5);   
    hess_app=zeros(5,5);
    for i=1:length(xi)        
        [tmp_grad, tmp_hess, tmp_hess_app]=get_grad_hessian(xi(i),yi(i),zi(i),ti(i),M);    

        grad = grad + tmp_grad;
        hess = hess + tmp_hess;   
        hess_app = hess_app + tmp_hess_app;
    end        
    err(step) = norm(predict(xi,yi,zi,M)-ti);    
    err2(step) = norm(predict(xi,yi,zi,M2)-ti);
    if(step > 1 && err(step) > err(step-1))
        break;
    end
    if(use_app)
        M= M- pinv(hess_app)*grad;
        M2= M2-pinv(hess)*grad;
    else
        M= M- pinv(hess)*grad;
        M2= M2-pinv(hess_app)*grad;
    end
end
% plot(xi,yi,'b.'); hold on;
% plot(M(2),M(3),'ro'); hold off;
% 1;

% function test_hessian()
% M=[0, 0, 0, 0, 3]';
% xi=100; yi=100; zi=0;  ti=50;
% [grad, hess] = get_grad_hessian(xi,yi,zi,ti,M);
% err0 = (predict(xi,yi,zi,M)-ti)^2;
% dM=[0,0,1,0,1]'*1E-2;
% err1 = (predict(xi,yi,zi,M+dM)-ti)^2;
% [err1-err0, grad'*dM, grad'*dM + 0.5*dM'*hess*dM];


function [grad,hess,hess_approx] = get_grad_hessian(xi,yi,zi,ti,M)
ts=M(1); xs=M(2); ys=M(3); zs=M(4); v=M(5);
R=sqrt((xs-xi)^2 + (ys-yi)^2 + (zs-zi)^2);
nx=(xs-xi)/R;
ny=(ys-yi)/R;
nz=(zs-zi)/R;
e=(ts+R/v-ti);

grad = e*[1 nx/v ny/v nz/v -R/v^2]';
sen = [1 nx/v ny/v nz/v -R/v^2]';

hess=[1 nx/v ny/v nz/v -R/v^2
    0  1/v*(1/v*nx^2+e/R*(1-nx^2))  1/v*(1/v*nx*ny-e/R*nx*ny) 1/v*(1/v*nx*nz-e/R*nx*nz) -1*nx/v^2*(R/v+e)
    0  0  1/v*(1/v*ny^2+e/R*(1-ny^2))  1/v*(1/v*ny*nz-e/R*ny*nz) -1*ny/v^2*(R/v+e)
    0  0   0  1/v*(1/v*nz^2+e/R*(1-nz^2))  -1*nz/v^2*(R/v+e)
    0 0 0 0 1*R/v^3*(R/v+2*e)];

hess = hess + hess';
for k=1:5
    hess(k,k)=hess(k,k)/2;
end
hess_approx = sen*sen'; % approximate hess

function r=predict(xi,yi,zi,M)

ts=M(1); xs=M(2); ys=M(3); zs=M(4); v=M(5);
r=xi*0;
for i=1:length(xi)
    R=sqrt((xs-xi(i))^2 + (ys-yi(i))^2 + (zs-zi(i))^2);
    r(i)=(ts+R/v);
end
    

