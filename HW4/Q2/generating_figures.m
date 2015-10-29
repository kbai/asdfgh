ii = 1:0.1:70;
mu = 35;
sig = 6;
yy = 1/(sqrt(2*pi)*sig)*exp(-(ii-mu).^2/(2*sig^2));
figure(1)
plot(ii,yy)
xlabel('p value');
ylabel('probability density');
title('Prior from geologist');
print('Geo_prior.pdf','-dpdf');


figure(2)
yy2 = 1./(ii);
plot(ii,yy2)
xlabel('p value');
ylabel('probability density');
title('Prior from physicist');
print('Phy_prior.pdf','-dpdf');


figure(3)
yy2yy = yy2.*yy;
plot(ii,yy2yy)
xlabel('p value');
ylabel('probability density');
title('Multiply two priors');
print('MUl_prior.pdf','-dpdf');
