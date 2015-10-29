xmin = 1;
xmax = 2;
delta = 0.001;
sum = 0;
for x = xmin:delta:xmax-delta
    sum = sum + 0.5*(f(x) + f(x+delta))*delta;
end