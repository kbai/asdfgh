function hw3_p1_d2

d = 5;
L = -1:0.1:11;

p = 15/2/sqrt(2*pi) .* exp(-50*(L-d).^2) + 3/2/sqrt(2*pi).* exp(-50*(L-d+1).^2);

figure(2)
plot(L,p,'-');
set(gca,'XTick',-1:2:11);    % get rid of the ticks
set(gca,'XTickLabel',str2mat('d_i-6','d_i-4','d_i-2','d_i','d_i+2','d_i+4','d_i+6'));
xlabel('L (cm)');
ylabel('P(L|d_i)');