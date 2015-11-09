%%%problem 1d
xi = [0 11 15 6 -7 3]';
yi = [0 0 6 13 10 -7]';
ui = [0.103 0.162  0.065  0.036 0.025 0.169]';
M0 = [8.1371 -5.1421 11.5066 30.3461]';  %initial guess 
X=zeros(1000,1);
P=zeros(1000,1);
R=zeros(1000,1);
for index=1:1:1
Minit = M0+1.0*[randn();0;0;randn()];
[con1,con2,Ms] = nonlinear_solver(xi,yi,ui,Minit); 
X(index) = Minit(1);
P(index) = Minit(4);
R(index) = con2;
end




