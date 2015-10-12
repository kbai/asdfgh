% fit a line to the data points, minimize the error in LS sense

function fit_data(file)

% default file name
if nargin < 1
    file = 'data.txt'; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load data

data = load(file); 
x = data(:,1);
y = data(:,2);
y0 = data(:,3);
% Q: how to read "data_h.txt"?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit a line
[k,b] = fit_a_line(x,y);

%
y_pred = k * x + b;

% Q: how to calculate the error of your fitting?

% plot
figure(1);
plot(x,y,'bo');
hold on;
plot(x,y_pred,'r*-');
plot(x,y0,'k-');
legend('data','pred','true','location','northwest');
hold off;
set(gca,'fontsize',14,'fontname','times');
xlabel('X');
ylabel('Y');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit a line, another method
[k2,b2] = fit_a_line2(x,y);

% output the results of two methods to screen
disp('k        b');
disp([k b; k2 b2]);

% plot 2d data
figure(2)
imagesc([k b; k2 b2]);
colorbar('vert');

end

% Given x, y data, fit a line
% use matlab function polyfit
function [k,b] = fit_a_line(x,y)
p = polyfit(x,y,1); 
k = p(1);
b = p(2);
end

% Given x,y data, fit a line
% use least square solution to Ax = b
function [k,b] = fit_a_line2(x,y)
N = length(x);
A = [x(:) ones(N,1)];

y = y(:);
p = (A'*A)\(A'*y);
k = p(1);
b = p(2);
end
