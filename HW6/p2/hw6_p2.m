function hw6_p2
set(0,'defaulttextfontname','times','defaulttextfontsize',14);
set(0,'defaultaxesfontname','times','defaultaxesfontsize',14);


% x{1} = [1 -3 4 5]';
% x{2} = [-0.1 0.3 -0.4 0.5]';
% x{3} = [101 97 104 105]';

x{1} = [1 2 3 4]';
x{2} = x{1} * 0.1;
x{3} = [10 12 14 16]';
%x{2} = x{1} * 0.1;
%x{3} = x{1} * 10;

m1 = 1; m2 = 1;
for i = 1:3
    y{i} = m1 + m2.*x{i} + [0;-0.15;0.1;-0.2];
end


figure(1);
for i = 1:3
    subplot(1,3,i)
    plot(x{i},y{i},'o');
    hold on;
    plot(x{i},m1 + m2.*x{i},'k-');
    hold all;
end

%
for i = 1:3
    G{i} = [ones(4,1) x{i}];
end

for i = 1:3
    disp(i);
    [U,SING,V] = svd(G{i})
    [EIGVEC,EIGV] = eig(G{i}'*G{i})
    
end


% test, add same level of error


figure(1)

for i = 1:3
   m = (G{i}'*G{i})^-1*G{i}'*y{i}
end
