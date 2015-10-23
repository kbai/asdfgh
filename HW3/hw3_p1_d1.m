function hw3_p1_d1
L = 5;
d = 3:0.01:7;

p = 15/2/sqrt(2*pi) .* exp(-50*(d-L).^2) + 5/2/sqrt(2*pi).* exp(-50*(d-L-1).^2);

figure(1)
plot(d,p,'-');
set(gca,'XTick',3:1:7);    % get rid of the ticks
set(gca,'XTickLabel',str2mat('L-2','L-1','L','L+1','L+2'));
xlabel('d_i (cm)');
ylabel('P(d_i|L)');
