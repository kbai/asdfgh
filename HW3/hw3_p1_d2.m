function hw3_p1_d2

d = 5;
L = 3:0.01:7;

p = 15/2/sqrt(2*pi) .* exp(-50*(L-d).^2) + 5/2/sqrt(2*pi).* exp(-50*(L-d+1).^2);

figure(2)
plot(L,p,'-');
set(gca,'XTick',3:1:7);    % get rid of the ticks
set(gca,'XTickLabel',str2mat('d_i-2','d_i-1','d_i','d_i+1','d_i+2'));
xlabel('L (cm)');
ylabel('P(L|d_i)');