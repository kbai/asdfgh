function hw3_p1_d3

L = 6:0.01:10;

f1 = 3*exp(-50*(9.1-L).^2) + exp(-50*(8.1-L).^2);
f2 = 3*exp(-50*(L-8.3).^2) + exp(-50*(L-7.3).^2);
f = f1 .* f2;
plot(L,f,'k-');

[~,loc] = max(f);
hold on;
plot(L(loc),f(loc),'r*');
text(L(loc),f(loc),num2str(L(loc),'%5.3f'));
hold off;

xlabel('L (cm)');
ylabel('const * P(L|d_1,d_2)');
end