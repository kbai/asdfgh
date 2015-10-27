xi = [0, 1, 2, 4]';
ti = [0 ,2, 3 ,10]';
M = [ 100, 100, 100]';
for ii = 1:1:10
    M = compute_gradient_hess(xi,ti,M);
end