function hw2_2
set(0,'defaulttextfontname','times','defaulttextfontsize',14);
set(0,'defaultaxesfontname','times','defaultaxesfontsize',14);

% load data
load ge118_hw2.mat
sigma = 0.1;   % data error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% least squares
disp('*** start least squares ***');
[m1_ls,m2_ls,mcov] = least_square(x,y,sigma);

disp('model covariance matrix:');
disp(mcov);
disp('standard deviation of m1 and m2:');
disp(sqrt(diag(mcov)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid search
disp('*** start monte carlo ***');
%m1 = -200:0.1:50;
%m2 = 1:0.1:10;
m1 = (m1_ls - 0.2):0.0005:(m1_ls + 0.2);
m2 = (m2_ls - 0.01):0.0001:(m2_ls + 0.01);

% L2 norm
[m1_l2,m2_l2,err_l2,err_all_l2] = grid_search(x,y,m1,m2,2);

% joint pdf
dm1 = m1(2) - m1(1);
dm2 = m2(2) - m2(1);
p_m1m2 = exp( -(err_all_l2-min(err_all_l2(:)))/(2*sigma^2) );
p_m1m2 = p_m1m2/(sum(p_m1m2(:))*dm1*dm2);

% marginal
p_m1 = dm2 * sum(p_m1m2,2);
p_m2 = dm1 * sum(p_m1m2,1);

m1 = m1(:);
m2 = m2(:);
p_m1 = p_m1(:);
p_m2 = p_m2(:);

figure;
subplot(121);
plot(m1,p_m1);
xlabel('m1'); ylabel('P(m1|d)');
subplot(122)
plot(m2,p_m2);
xlabel('m2'); ylabel('P(m2|d)');

% \sigma m1 m2
sigma_m1 = sqrt( dm1 * sum(p_m1 .* m1.^2) - (dm1 * sum(p_m1 .* m1)).^2 );
sigma_m2 = sqrt( dm2 * sum(p_m2 .* m2.^2) - (dm2 * sum(p_m2 .* m2)).^2 );
disp('standard deviation of m1 and m2:');
disp([sigma_m1 sigma_m2]);


end

% LEAST SQUARE SOLUTION
function [m1, m2, mcov] = least_square(xdata,ydata,sigma)
G = [ones(length(xdata),1) xdata(:)];
tmp = (G'*G)^(-1) * G' * ydata;
m1 = tmp(1);
m2 = tmp(2);
mcov = inv(G'*G/sigma^2);
end

% GRID SEARCH
function [m1_best,m2_best,err_best,err] = grid_search(xdata,ydata,m1,m2,flag)
% error over the model space
err = zeros(length(m1),length(m2));
for i = 1:length(m1)
    for j = 1:length(m2)
        err(i,j) = misfit(xdata,ydata,m1(i),m2(j),flag);
    end
end

% find the optimum solution
[err_best,id] = min(err(:));
[I,J] = ind2sub(size(err),id);
m1_best = m1(I);
m2_best = m2(J);
end

% CALCULATE THE MISFIT
function err = misfit(xdata,ydata,m1,m2,flag)
ypred = m1 + m2 * xdata;
if flag == 1       % L1 norm
    err = sum(abs(ypred-ydata));
elseif flag == 2   % L2 norm
    err = sum((ypred-ydata).^2);
end
end