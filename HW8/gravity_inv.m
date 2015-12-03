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

% plot the data
figure(11)
subplot(211)
plot(xd, di, 'ko-');
xlabel('x (m)'); ylabel('\muGal'); xlim([-5 15]);
title('data');
subplot(212)
plot(xi, LS_Dm);
xlabel('x (m)'); ylabel(' \Deltam (kg)'); xlim([-5 15]);
title('standard least square');
print('figure/e_result.png','-dpng');

% SVD of GG
[UU,DD,VV]=svd(GG);

figure(2);
plot(1:1:min(length(UU),length(VV)),log(diag(DD))/log(10));
ylabel('log_1_0 (S)');
title('singular values');
print('figure/e_snv.png','-dpng');

Re = eye(Lm);

%% part (f)

[GI_Dm,Rf] = generalized_inverse(di,GG);   %generalized inverse without truncation

figure(3)

plot(xi,LS_Dm,'k-');
hold on
plot(xi,GI_Dm,'r-');
legend('LS','SVD');
print('figure/f.png','-dpng');

%% part (g)
p=10;
figure(4)

%for p = 10:5:20
 for p = 5:10   
[GI_Dm_trunc,Rg{p}] = generalized_inverse(di,GG,p);   %generalized inverse with truncation

plot(xi,GI_Dm_trunc);


hold on

end

%% part (h)


nm=[];
ne=[];
alpha_all = 10.^[-10:-1];
for i= 1:length(alpha_all);
alpha_s = alpha_all(i);
Tik_Dm = (GG'*GG + alpha_s^2 * eye(Lm))\(GG'*di);
nm = [nm,norm(Tik_Dm)];
ne = [ne,norm(GG*Tik_Dm-di)];

figure(5)
plot(xi,Tik_Dm);
hold on

end


figure(6)
loglog(ne,nm,'o-');
loglog(ne(5),nm(5));

xlabel('Residual norm |Gm-d|_2');
ylabel('Solution norm |m|_2')

%% part(i)

L1 = diag(ones(1,Lm),0) + diag(-1*ones(1,Lm-1),1);
L1(end,:) = [];


alpha = 1e-5;

Tik1_Dm = (GG'*GG + alpha^2 * (L1'*L1))\(GG'*di);

Ri1 = (GG'*GG + alpha^2 * (L1'*L1))\(GG'*GG);

figure(7)

plot(xi,Tik1_Dm);

L2 = diag(ones(1,Lm),0) + diag(-2*ones(1,Lm-1),1) + diag(ones(1,Lm-2),2);
L2(end,:) = [];
L2(end,:) = [];

Tik2_Dm = (GG'*GG + alpha^2 * (L2'*L2))\(GG'*di);

Ri2 = (GG'*GG + alpha^2 * (L2'*L2))\(GG'*GG);

hold on

plot(xi,Tik2_Dm);

%% part(j)





    
        












