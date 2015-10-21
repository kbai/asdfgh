function hw2p1()
rng(0); % this generate determined random variable
xi=[0 10 15 6 -7 3]';
yi=[0 0 6 13 10 7]';
zi=[0 0 0 0 0 0]';
ti=[322.418 321.031 321.228 323.093 324.415 322.706]';

M0=do_one_ti(xi,yi,zi,ti);
ti_0=predict(xi,yi,zi,M0);


for it=1:1000
    ti_new = ti + 0.01*randn(6,1);
    % inexact Hessian
    M(:,it)=do_one_ti(xi,yi,zi,ti_new);
end
scatter(M(2,:), M(3,:), 10, M(1,:)); colorbar;
hold on;plot(M0(2),M0(3),'rp','markersize',20,'markerfacecolor','r');
print -depsc fig01.eps
plot([0 10 15 6 -7 3],[0 0 6 13 10 7],'kv');
hold off;
print -depsc fig02.eps

function M=do_one_ti(xi,yi,zi,ti)

M=[300 20 -10 10 2]'; % initial guess
for step = 1:1000
    grad=zeros(5,1);
    hess=zeros(5,5);   
    for i=1:length(xi)
        [tmp_grad, tmp_hess]=get_grad_hessian(xi(i),yi(i),zi(i),ti(i),M);    
        grad = grad + tmp_grad;
        hess = hess + tmp_hess;   
    end        
    err(step) = norm(predict(xi,yi,zi,M)-ti);    
    if(step > 1 && err(step) > err(step-1))
        break;
    end
	M = M- inv(hess)*grad;
end
plot(err);

function [grad,hess] = get_grad_hessian(xi,yi,zi,ti,M)
ts=M(1); xs=M(2); ys=M(3); zs=M(4); v=M(5);
R=sqrt((xs-xi)^2 + (ys-yi)^2 + (zs-zi)^2);
nx=(xs-xi)/R;
ny=(ys-yi)/R;
nz=(zs-zi)/R;
e=(ts+R/v-ti);

grad = e*[1 nx/v ny/v nz/v -R/v^2]';
sen = [1 nx/v ny/v nz/v -R/v^2]';
hess = sen*sen'; % approximate hess

function r=predict(xi,yi,zi,M)

ts=M(1); xs=M(2); ys=M(3); zs=M(4); v=M(5);
r=xi*0;
for i=1:length(xi)
    R=sqrt((xs-xi(i))^2 + (ys-yi(i))^2 + (zs-zi(i))^2);
    r(i)=(ts+R/v);
end
