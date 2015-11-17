%calculate generalized inverses and svds
function hw7p2()
%problem 2
%To truncate these SVDs we will remove the smallest singular value (the one
%most likely responsible for propagation of large errors).

clc
close all

G1=[1 1 1 1; 1 -3 4 5]';
G2=[1 1 1 1; -0.1 0.3 -0.4 0.5]';
G3=[1 1 1 1; 101 97 104 105]';
G4=[1 1 1 1; 3 3 3 3]';

[U,S,V] = svd(G1);
[U,S,V] = svd(G2);
[U,S,V] = svd(G3);

R1=pseduoInverse(G1);
R2=pseduoInverse(G2);
R3=pseduoInverse(G3);
R4=pseduoInverse(G4);

T3 = GenInvTruncated(G3)
T1 = GenInvTruncated(G1)

%%%%%%%
%Part D
d1 = [10 11 9 12]';
d2 = [10.1 11.4 8.7 9.8]';

x1 = [1 -3 4 5]';
x2 = [-0.1 0.3 -0.4 0.5]';
x3 = [101 97 104 105]';
x4 = [3 3 3 3]';

G1x = linspace(-5,7);
G3x = linspace(95,107);

%calculate generalized inverse model solution
%using ORIGINAL SVD
m1d1 = R1*d1
m1d2 = R1*d2
m3d1 = R3*d1
m3d2 = R3*d2

%using TRUNCATED SVD
mt1d1 = T1*d1
mt1d2 = T1*d2
mt3d1 = T3*d1
mt3d2 = T3*d2

figure(1)                   %G1 and d1
plot(x1,d1,'k^')
plotOriginal(m1d1,G1x)
plotTruncated(mt1d1,G1x)
title('G1 - data1')
xlabel('x')
ylabel('d')
legend('data','original','truncated','Location','SouthEast')

figure(2)                   %G1 and d2
plot(x1,d2,'k^')
plotOriginal(m1d2,G1x)
plotTruncated(mt1d2,G1x)
title('G1 - data2')
xlabel('x')
ylabel('d')
legend('data','original','truncated','Location','SouthEast')

figure(3)                   %G3 and d1
plot(x3,d1,'k^')
plotOriginal(m3d1,G3x)
plotTruncated(mt3d1,G3x)
title('G3 - data1')
xlabel('x')
ylabel('d')
legend('data','original','truncated')

figure(4)                   %G3 and d2
plot(x3,d2,'k^')
plotOriginal(m3d2,G3x)
plotTruncated(mt3d2,G3x)
title('G3 - data2')
xlabel('x')
ylabel('d')
legend('data','original','truncated')

function plotOriginal(m,Gx)
hold on
y = m(1) + m(2)*Gx;
plot (Gx,y,'r')


function plotTruncated(m,Gx)
hold on
y = m(1) + m(2)*Gx;
plot (Gx,y,'b')


%GENERALIZED INVERSE FUNCTION
function r=pseduoInverse(G)
[U,S,V]=svd(G); % check if U*S*V' == G
eps=1E-16;
for i=1:min(size(S))
    if(S(i,i) < eps)
        S(i,i) = 0;
    else
        S(i,i) = 1.0/S(i,i);
    end
end
r=V*S'*U';

%TRUNCATE FUNCTION
function r=GenInvTruncated(G)
[U,S,V]=svd(G);     %find SVD
eps=1E-16;
for i=1:min(size(S))
    if(S(i,i) < eps)
        S(i,i) = 0;
    else
        S(i,i) = 1.0/S(i,i);
    end
end
%now truncate the smallest singular value (in this case it's S3(2,2))
Struncated = S;
Struncated(2,2) = 0;    %in this case each of the G matrices only have 2 singular values
%recalculate the generalized inverse
r = V * Struncated' * U';





