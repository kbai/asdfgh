function hw3_p1_d1
L = 5;
d = -1:0.1:11;

p = 15/2/sqrt(2*pi) .* exp(-50*(d-L).^2) + 3/2/sqrt(2*pi).* exp(-50*(d-L-1).^2);

figure(1)
plot(d,p,'-');
set(gca,'XTick',-1:2:11);    % get rid of the ticks
set(gca,'XTickLabel',str2mat('L-6','L-4','L-2','L','L+2','L+4','L+6'));
xlabel('d_i (cm)');
ylabel('P(d_i|L)');
