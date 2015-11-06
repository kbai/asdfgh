%%%problem 1d
xi = [0 11 15 6 -7 3]';
yi = [0 0 6 13 10 -7]';
ui = [0.103 0.162  0.065  0.036 0.025 0.169]';
M0 = [8 -5 10 30]';  %initial guess 
M0 = M0+1.0*randn(4,1);
Ms = nonlinear_solver(xi,yi,ui,M0); 




