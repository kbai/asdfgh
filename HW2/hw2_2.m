function hw2_2
set(0,'defaulttextfontname','times','defaulttextfontsize',14);
set(0,'defaultaxesfontname','times','defaultaxesfontsize',14);

% load data
load ge118_hw2.mat

% plot data
figure(1)
plot(x,y,'k*'); 
grid on; axis equal;
xlabel('x'); ylabel('y');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid search
m1 = -150:0.1:0;
m2 = 1:0.01:10;
% L1 norm
[m1_l1,m2_l1,err_l1,err_all_l1] = grid_search(x,y,m1,m2,1);
% L2 norm
[m1_l2,m2_l2,err_l2,err_all_l2] = grid_search(x,y,m1,m2,2);

% plotting
figure(2)
subplot(121)
pcolor(m2,m1,err_all_l1); 
hold on;
plot(m2_l1,m1_l1,'ko');
shading flat;
caxis(crange(err_all_l1));
xlabel('m2');ylabel('m1');title('L1 norm');
hold off; 

subplot(122)
pcolor(m2,m1,err_all_l2);
hold on;
plot(m2_l2,m1_l2,'ko');
shading flat;
caxis(crange(err_all_l2));
xlabel('m2');ylabel('m1');title('L2 norm');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% least square
[m1_ls,m2_ls] = least_square(x,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('L1 norm           L2 norm             LS');
disp('m1');
disp([m1_l1  m1_l2  m1_ls]);
disp('m2');
disp([m2_l1  m2_l2  m2_ls]);

end

% Least square solution
% input: x y data
% output: optimum parameters m1 m2
function [m1, m2] = least_square(x,y)
G = [ones(length(x),1) x(:)];
tmp = (G'*G)^(-1) * G' * y;
m1 = tmp(1);
m2 = tmp(2);
end

% Grid search
% input:
%    x y   - data; 
%    m1 m2 - range for search; 
%    flag  - norm for misfit
% output: 
%     m1_best m2_best err_best - optimum parameters and its misift
%     err - misfit of each point in model space N(m1) x N(m2)
function [m1_best,m2_best,err_best,err] = grid_search(x,y,m1,m2,flag)
err = zeros(length(m1),length(m2));
for i = 1:length(m1)
    for j = 1:length(m2)
        y_pred = m1(i) + m2(j) * x;
        err(i,j) = misfit(y,y_pred,flag);
    end
end

[err_best,id] = min(err(:));
[I,J] = ind2sub(size(err),id);
m1_best = m1(I);
m2_best = m2(J);
end

% calculate the misfit
function err = misfit(y0,y1,flag)
if flag == 1       % L1 norm
    err = sum(abs(y0-y1));
elseif flag == 2   % L2 norm
    err = sum((y0-y1).^2);
end
end

% caxis for error plot
function vec = crange(err)
minval = min(err(:));
maxval = minval + 0.15 * (max(err(:)) - minval);
vec = [minval maxval];
end