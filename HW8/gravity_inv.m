
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
xi = [0:1:100]*0.1;

unit = 1e-8; %unit of microgals 1e-8m/s^2

GG = (1/unit)*G*D./(D^2+(X'*ones(1,Lm) - ones(Ld,1)*xi).^2).^1.5;

Dgz = GG*m;

figure(1);
plot(X,Dgz);

%% part (e)
load('ge118_hw8.mat');
Ld = length(xd);

Lm = length(m);

GG = (1/unit)*G*D./(D^2+(xd*ones(1,Lm) - ones(Ld,1)*xi).^2).^1.5;

[UU,DD,VV]=svd(GG);

figure(2);

plot(1:1:min(length(UU),length(VV)),log(diag(DD))/log(10));

LS_Dm = (GG'*GG)\(GG'*di);

Re = eye(Lm);

%% part (f)

[GI_Dm,Rf] = generalized_inverse(di,GG);   %generalized inverse without truncation

figure(3)

plot(1:1:length(Lm),LS_Dm);

hold on

plot(1:1:length(Lm),GI_Dm);

%% part (g)
p=10;
figure(4)

for p = 10:5:20
    
[GI_Dm_trunc,Rg{p}] = generalized_inverse(di,GG,p);   %generalized inverse with truncation

plot(xi,GI_Dm_trunc);


hold on

end

%% part (h)

figure(5)
nm=[];
ne=[]
for alpha_s =[1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12]

Tik_Dm = (GG'*GG + alpha_s * eye(Lm))\(GG'*di);
nm = [nm,norm(Tik_Dm)];
ne = [ne,norm(GG*Tik_Dm-di)];

plot(xi,Tik_Dm);

hold on

end
figure(6)

loglog(ne,nm,'o-');

xlabel('Residual norm |Gm-d|_2');
ylabel('Solution norm |m|_2')

%% part(i)

L1 = zeros(Lm-1,Lm);

for ii = 1:1:Lm-1
    L1(ii,ii) = -1;
    L1(ii,ii+1) = 1;
end

alpha = 1e-10;

Tik1_Dm = (GG'*GG + alpha * (L1'*L1))\(GG'*di);

Ri1 = (GG'*GG + alpha * (L1'*L1))\(GG'*GG);

figure(7)

plot(xi,Tik1_Dm);

L2 = zeros(Lm-2,Lm);

for ii = 1:1:Lm-2
    L2(ii,ii) = -1;
    L2(ii,ii+1) = 2;
    L2(ii,ii+2) = -1;
end

Tik2_Dm = (GG'*GG + alpha * (L2'*L2))\(GG'*di);

Ri2 = (GG'*GG + alpha * (L2'*L2))\(GG'*GG);

hold on

plot(xi,Tik2_Dm);

%% part(j)





    
        












