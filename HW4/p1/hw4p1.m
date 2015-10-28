function hw4p1
set(0,'defaulttextfontname','times','defaulttextfontsize',14);
set(0,'defaultaxesfontname','times','defaultaxesfontsize',14);

sigma = 0.001;
std_mc = [0.099670 0.137469 0.149719 0.691071]';
%%%hw 2 problem 1d
xi = [0 11 15 6 -7 3]';
yi = [0 0 6 13 10 -7]';
ui = [0.103 0.162  0.065  0.036 0.025 0.169]';
M0 = [8 -5 10 30]';  %initial guess, xs ys zs P
[Ms, mcov] = nonlinear_solver(xi,yi,ui,M0,sigma); 

% model covariance matrix
disp('model covariance matrix:');
disp(mcov);
disp('standard deviation');
std_m = sqrt(diag(mcov))

% calculate correlation matrix
s = diag(std_m);
disp('coefficient matrix:');
mcoef = s^(-1)*mcov*s^(-1)

return;
%%%hw 2 problem 1e
Mrec = zeros(4,1000);

for it = 1:1:1000
    uerror = ui + 0.001*randn(6,1);
    Merror = nonlinear_solver(xi,yi,uerror,M0);
    Mrec(:,it) = Merror;
end

colormap(jet);
scatter(Mrec(1,:),Mrec(2,:),30, Mrec(3,:),'fill');
stdx=std(Mrec(1,:));
stdy=std(Mrec(2,:));
stdz=std(Mrec(3,:));
stdp=std(Mrec(4,:));
c = colorbar;
ylabel(c,'depth(km)');
xlabel('x(km)');
ylabel('y(km)');

print('measurements_error.pdf','-dpdf');
end
