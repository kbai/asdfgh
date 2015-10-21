function hw3()
xi=[0 10 15 6 -7 3]';
yi=[0 0 6 13 10 7]';
zi=[0 0 0 0 0 0]';
ti=[322.418 321.031 321.228 323.093 324.415 322.706]';

sigma0=zeros(6,1)+0.01;
M0=do_one_ti(xi,yi,zi,ti, sigma0);
ti_0=predict(xi,yi,zi,M0);

sigma1=[0.007E-3, 0.01, 0.01, 0.007E-3, 0.1, 0.01]';
M1=do_one_ti(xi,yi,zi,ti, sigma1);
ti_1=predict(xi,yi,zi,M1);

sigma2=[0.007, 90000.01, 0.01, 90000.007, 90000.1, 0.01]';
M2=do_one_ti(xi,yi,zi,ti, sigma2);
ti_2=predict(xi,yi,zi,M2);

plot(xi,yi,'bv'); hold on;
plot(xi([2,4,5]),yi([2,4,5]),'r.');
plot(M0(2),M0(3),'kp','markersize',14,'markerfacecolor','k');
plot(M1(2),M1(3),'bp','markersize',10,'markerfacecolor','b');
plot(M2(2),M2(3),'rp','markersize',10,'markerfacecolor','r');
hold off;
print -depsc2 fig1.eps
[M0, M1, M2]
[ti_0-ti, ti_1-ti, ti_2-ti]


function M=do_one_ti(xi,yi,zi,ti,sigma)

M=[300 20 -10 10 2]'; % initial guess
for step = 1:1000
    grad=zeros(5,1);
    hess=zeros(5,5);   
    for i=1:length(xi)
        [tmp_grad, tmp_hess]=get_grad_hessian(xi(i),yi(i),zi(i),ti(i),M);    
        grad = grad + tmp_grad/sigma(i)^2;
        hess = hess + tmp_hess/sigma(i)^2;   
    end        
    err(step) = norm(predict(xi,yi,zi,M)-ti);    
    if(step > 1 && err(step) > err(step-1))
        break;
    end
	M = M- pinv(hess)*grad;
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
