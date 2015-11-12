function hw6_p2
set(0,'defaulttextfontname','times','defaulttextfontsize',14);
set(0,'defaultaxesfontname','times','defaultaxesfontsize',14);


x1 = [1 -3 4 5];
x2 = [-0.1 0.3 -0.4 0.5];
x3 = [101 97 104 105];


m1 = 1; m2 = 1;

y1 = m1 + m2.*x1 + randn(1,4);
y2 = m1 + m2.*x2 + randn(1,4);
y3 = m1 + m2.*x3 + randn(1,4);

figure(1);
subplot(131)
plot(x1,y1,'ro'); grid on; xlabel('x');ylabel('y'); axis square;
subplot(132)
plot(x2,y2,'go'); grid on; xlabel('x');ylabel('y'); axis square;
subplot(133)
plot(x3,y3,'bo'); grid on; xlabel('x');ylabel('y'); axis square;

%
G{1} = [ones(4,1) x1(:)];
G{2} = [ones(4,1) x2(:)];
G{3} = [ones(4,1) x3(:)];

for i = 1:3
    disp(i);
    [U,SING,V] = svd(G{i})
    [EIGVEC,EIGV] = eig(G{i}'*G{i})
    
end
