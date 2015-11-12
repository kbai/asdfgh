function hw6p2()
xi=[0 10 15 6 -7 3]';
yi=[0 0 6 13 10 7]';
zi=[0 0 0 0 0 0]';
ti=[322.418 321.031 321.228 323.093 324.415 322.706]';
M0=[300 10 -10 10 5]';

lambda=zeros(1000,1);
[M1,err1]=do_one_ti(xi,yi,zi,ti, M0,lambda);
lambda(1:3)=10;
[M2,err2]=do_one_ti(xi,yi,zi,ti, M0,lambda);

% the error contour is ploted assuming ys,zs,v is at the best guess
Mbest = M1(:,end);
perturb_ts = linspace(-30,30,100);
perturb_xs = linspace(-30,30,200);
for k=1:length(perturb_ts)
    for j=1:length(perturb_xs)
        Mnew=Mbest;
        Mnew(1) = Mnew(1) + perturb_ts(k);
        Mnew(2) = Mnew(2) + perturb_xs(j);
        F(k,j) = 0.5*norm(predict(xi,yi,zi,Mnew)-ti)^2;
    end
end

% be careful about x and y axis
contour(Mbest(1) + perturb_ts, Mbest(2)+perturb_xs, F');
hold on;
plot(M1(1,:), M1(2,:),'ro-','linewidth',2,'markersize',7); 
plot(M2(1,:), M2(2,:),'bo-','linewidth',2,'markersize',7); 
axis equal
legend({'','GN','LM'})
hold off;
print -depsc LM.eps





function [Msave,err]=do_one_ti(xi,yi,zi,ti, M0, lambda)

M=M0;
Msave=M;
for step = 1:10
    grad=zeros(5,1);
    hess=zeros(5,5);     % exact hessian
    hess_app=zeros(5,5); % approximate hessian
    for i=1:length(xi)        
        [tmp_grad, tmp_hess_app]=get_grad_hessian(xi(i),yi(i),zi(i),ti(i),M);
        grad = grad + tmp_grad;
        hess_app = hess_app + tmp_hess_app;
    end        
    err(step) = norm(predict(xi,yi,zi,M)-ti)/sqrt(length(xi)); 
    M = M - inv(hess_app + lambda(step)*diag(diag(hess_app)) ) *grad;
    
    Msave=[Msave,M];
end

function [grad,hess_approx] = get_grad_hessian(xi,yi,zi,ti,M)
ts=M(1); xs=M(2); ys=M(3); zs=M(4); v=M(5);
R=sqrt((xs-xi)^2 + (ys-yi)^2 + (zs-zi)^2);
nx=(xs-xi)/R;
ny=(ys-yi)/R;
nz=(zs-zi)/R;
e=(ts+R/v-ti);

grad = e*[1 nx/v ny/v nz/v -R/v^2]';
sen = [1 nx/v ny/v nz/v -R/v^2]';
hess_approx = sen*sen'; % approximate hess

function r=predict(xi,yi,zi,M)
ts=M(1); xs=M(2); ys=M(3); zs=M(4); v=M(5);
r=xi*0;
for i=1:length(xi)
    R=sqrt((xs-xi(i))^2 + (ys-yi(i))^2 + (zs-zi(i))^2);
    r(i)=(ts+R/v);
end
    