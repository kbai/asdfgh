function [ M_new ] = compute_gradient_hess( x,t,M )

A=M(1);
B=M(2);
C=M(3);

misfit = t - (A + B*x +C*x.^2);
G = [ ones(length(x),1), x ,x.^2];
H = G'*G;

M_new = M + pinv(H) * G'* misfit;

norm(misfit)
end

