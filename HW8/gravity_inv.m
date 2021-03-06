set(0,'defaulttextfontname','times','defaulttextfontsize',14);
set(0,'defaultaxesfontname','times','defaultaxesfontsize',14);
%% part (d)
G = 6.67e-11; %gravitational constant
m = zeros(101,1);
D = 2;

m(1) = 500;    % set mass anomaly at 1 , 31 and 91
m(31) = 300;
m(91) = 100;

X = -10:0.1:20;  %measurement location

Lm = length(m);
Ld = length(X);
xi = [0:1:100]*0.1;   % xi: location of the anomalies mi

unit = 1e-8; %unit of microgals 1e-8m/s^2

GG = (1/unit)*G*D./(D^2+(X'*ones(1,Lm) - ones(Ld,1)*xi).^2).^1.5;

Dgz = GG*m;

figure(1);
subplot(211)
plot(X,Dgz); xlim([-10 20]);
xlabel('x (m)'); ylabel('\muGal');
hold on;
subplot(212);
plot(xi,m); xlim([-10 20]);
xlabel('x (m)'); ylabel('\Deltamkg');
print('figure/d.png','-dpng');

%% part (e)
load('ge118_hw8.mat');
Ld = length(xd);

Lm = length(m);

GG = (1/unit)*G*D./(D^2+(xd*ones(1,Lm) - ones(Ld,1)*xi).^2).^1.5;

% Least square solution
LS_Dm = (GG'*GG)\(GG'*di);

%model prediction
d_pred = GG*LS_Dm;
% plot the data
figure(11)
subplot(211)
plot(xi, LS_Dm);
xlabel('x (m)'); ylabel(' \Deltam (kg)'); xlim([-5 15]);
title('standard least square');

subplot(212)
plot(xd, di, 'ko-');
hold on
plot(xd,d_pred,'b*-');
legend('data','model prediction');
xlabel('x (m)'); ylabel('\muGal'); xlim([-5 15]);
title('data');
print('figure/e_result.png','-dpng');

% SVD of GG
[UU,DD,VV]=svd(GG);

figure(2);
plot(1:1:min(length(UU),length(VV)),log(diag(DD))/log(10));
ylabel('log_1_0 (S)');
title('singular values');
print('figure/e_snv.png','-dpng');


Ginv_LLS = (GG'*GG)\GG';

Re_LLS = Ginv_LLS*GG;  %resolution matrix for least squares


%% part (f)

[GI_Dm,Rf] = generalized_inverse(di,GG);   %generalized inverse without truncation
d_predf = GG*GI_Dm;
figure(3)
subplot(211)

plot(xi,LS_Dm,'k-');
hold on
plot(xi,GI_Dm,'r-');
legend('LS','SVD');

subplot(212)
plot(xd, di, 'ko-');
hold on
plot(xd,d_predf,'b*-');
legend('data','model prediction');


print('figure/f.png','-dpng');

%% part (g)
nm=[];
ne=[];
p=10;
figure(4)
%for p = 10:5:20
for p = 5:30
    [GI_Dm_trunc,Rg{p}] = generalized_inverse(di,GG,p);   %generalized inverse with truncation
    plot(xi,GI_Dm_trunc);
    hold on
    ne = [ne,norm(di-GG*GI_Dm_trunc)];
    nm = [nm,norm(GI_Dm_trunc)];
end

figure(20);
loglog(nm,ne,'-o'); 
text(nm,ne,strread(num2str(5:1:30),'%s'));
legend('truncation number')
xlabel('||m||^2');
ylabel('||d-Gm||^2');
print('figure/g_trunc.png','-dpng');
hold on
xlim([3150,3450]);
ylim([0.1345,0.1395]);
xlabel('||m||^2');
ylabel('||d-Gm||^2');
print('figure/g_trunc_zoomin.png','-dpng');

%% part (h)
nm=[];
ne=[];
alpha_all = 10.^[-10:-1];
for i= 1:length(alpha_all);
    Tik_Dm = (GG'*GG + alpha_all(i)^2 * eye(Lm))\(GG'*di);
    nm = [nm,norm(Tik_Dm)];
    ne = [ne,norm(GG*Tik_Dm-di)];
    
    if i>=5 && i<= 7
        figure(5)
        if i == 5
            plot(xi,Tik_Dm,'k-');
        elseif i == 6
            plot(xi,Tik_Dm,'r-');
        elseif i == 7
            plot(xi,Tik_Dm,'b-');
        end
        hold on
        ylim([-1000 1000]);
    end
    Rt{i} = (GG'*GG + alpha_all(i)^2 * eye(Lm))\(GG'*GG);
end
xlabel('x (m)'); ylabel('\Deltamkg');
legend('\alpha=1e-5','\alpha=1e-4','\alpha=1e-3');
print('figure/h_solutions.png','-dpng');


figure(6)
loglog(ne,nm,'o-');
hold on;
plot(ne(6),nm(6),'r*');
xlabel('Residual norm |Gm-d|_2');
ylabel('Solution norm |m|_2');
hold off;
title('find alpha - 0 order');
print('figure/h_alphas.png','-dpng');

%% part(i)
% first order Tikhonov
L1 = diag(ones(1,Lm),0) + diag(-1*ones(1,Lm-1),1);
L1(end,:) = [];


% find the best alpha
alpha_all = 10.^[-10:-1];
nm = [];
ne = [];
for i = 1:length(alpha_all)
    Tik_Dm = (GG'*GG + alpha_all(i)^2 * (L1'*L1))\(GG'*di);
    nm = [nm,norm(L1*Tik_Dm)];
    ne = [ne,norm(GG*Tik_Dm-di)];
    Rt{i} = (GG'*GG + alpha_all(i)^2 * (L1'*L1))\(GG'*GG); 
end
figure(31);
subplot(121)
loglog(ne,nm,'o-');
hold on;
plot(ne(7),nm(7),'r*');
text(ne(7),nm(7),'\alpha=1e-4');
title('find alpha - first order');
xlabel('Residual norm |Gm-d|_2');
ylabel('Solution norm |Lm|_2');


alpha = 1e-4;
Tik1_Dm = (GG'*GG + alpha^2 * (L1'*L1))\(GG'*di);
Ri1 = (GG'*GG + alpha^2 * (L1'*L1))\(GG'*GG);   % resolution matrix



figure(7)
plot(xi,Tik1_Dm,'k-');

% second order Tikhonov

L2 = diag(ones(1,Lm),0) + diag(-2*ones(1,Lm-1),1) + diag(ones(1,Lm-2),2);
L2(end,:) = [];
L2(end,:) = [];

% find the best alpha
alpha_all = 10.^[-10:-1];
nm = [];
ne = [];
for i = 1:length(alpha_all)
    Tik_Dm = (GG'*GG + alpha_all(i)^2 * (L2'*L2))\(GG'*di);
    nm = [nm,norm(L2*Tik_Dm)];
    ne = [ne,norm(GG*Tik_Dm-di)];
end

figure(31);
subplot(122)
loglog(ne,nm,'o-');
hold on;
plot(ne(8),nm(8),'r*');
text(ne(8),nm(8),'\alpha=1e-3');
title('find alpha - second order');
xlabel('Residual norm |Gm-d|_2');
ylabel('Solution norm |Lm|_2');
print('figure/i_alphas.png','-dpng');

alpha = 1e-3;
Tik2_Dm = (GG'*GG + alpha^2 * (L2'*L2))\(GG'*di);

Ri2 = (GG'*GG + alpha^2 * (L2'*L2))\(GG'*GG);  % resolution matrix

figure(7)
hold on
plot(xi,Tik2_Dm,'r-');
xlabel('x (m)'); ylabel('\Deltamkg');
legend('first-order','second-order');
print('figure/i_solution.png','-dpng');


figure
plot(xd,di,'ko-');
hold on;
plot(xd,GG * Tik1_Dm,'r-');
plot(xd,GG * Tik2_Dm,'b-');
xlabel('x (m)'); ylabel('\muGal');
hold off;
legend('data','first order','second order');
print('figure/i_pred.png','-dpng');

%% part(j)
figure
imagesc(Re_LLS); colorbar('vert');
title('(e) Least square');
print('figure/j1.png','-dpng');

figure
imagesc(Rf); colorbar('vert');
title('(f) Generalized inverse');
print('figure/j2.png','-dpng');

figure
imagesc(Rg{11}); colorbar('vert');
title('(g) p = 11');
print('figure/j3.png','-dpng');

figure
imagesc(Rg{16}); colorbar('vert');
title('(g) p = 16');
print('figure/j4.png','-dpng');

figure
imagesc(Rg{21}); colorbar('vert');
title('(g) p = 21');
print('figure/j5.png','-dpng');

figure
imagesc(Rt{6}); colorbar('vert');
title('\alpha = 1e-5');
print('figure/j6.png','-dpng');

figure
imagesc(Rt{7}); colorbar('vert');
title('\alpha = 1e-4');
print('figure/j7.png','-dpng');

figure
imagesc(Rt{8}); colorbar('vert');
title('\alpha = 1e-3');
print('figure/j8.png','-dpng');

figure
imagesc(Ri1); colorbar('vert');
title('(i) first-order');
print('figure/j9.png','-dpng');

figure
imagesc(Ri2); colorbar('vert');
title('(i) second-order');
print('figure/j10.png','-dpng');





    
        












