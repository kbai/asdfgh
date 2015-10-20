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
xlim([-250 250]); ylim([-250 250]);
set(gca,'XTick',[-250:50:250]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% grid search
m1 = -200:0.1:50;
m2 = 1:0.01:10;

% L1 norm
[m1_l1,m2_l1,err_l1,err_all_l1] = grid_search(x,y,m1,m2,1);

% L2 norm
[m1_l2,m2_l2,err_l2,err_all_l2] = grid_search(x,y,m1,m2,2);

% plotting
figure(2)
colormap(jet);
subplot(121)
pcolor(m2,m1,err_all_l1);  % error map
hold on;
plot(m2_l1,m1_l1,'go');    % optimum solution
shading flat;
caxis(crange(err_all_l1));
colorbar('horiz');
xlabel('m2');ylabel('m1');title('L1 norm');
hold off; 

subplot(122)
pcolor(m2,m1,err_all_l2);  % error map
hold on;
plot(m2_l2,m1_l2,'go');    % optimum solution
shading flat;
caxis(crange(err_all_l2));
colorbar('horiz');
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

% LEAST SQUARE SOLUTION
function [m1, m2] = least_square(xdata,ydata)
G = [ones(length(xdata),1) xdata(:)];
tmp = (G'*G)^(-1) * G' * ydata;
m1 = tmp(1);
m2 = tmp(2);
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

% caxis for error plot
function vec = crange(err)
minval = min(err(:));
maxval = minval + 0.15 * (max(err(:)) - minval);
vec = [minval maxval];
end