%%%problem 2f
xi = [0 11 15 6 -7 3]';
yi = [0 0 6 13 10 -7]';
ui = [0.103 0.162  0.065  0.036 0.025 0.169]';
M0 = [10 -10 5 30]';  %initial guess 
Ms = nonlinear_solver(xi,yi,ui,M0); 





