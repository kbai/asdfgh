t = [0 1 2 4]';

G = [ones(4,1) t t.^2];
d = [0 2 3 10]';

m = (G'*G)^(-1)*G'*d;

plot(t,d,'ko');
hold on;
plot(t, m(1) + m(2) .* t + m(3) .* t.^2, 'r-');