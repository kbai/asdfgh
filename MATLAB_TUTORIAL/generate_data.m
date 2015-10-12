% Generate synthetic data points for LS fitting
% (1) generate a line; 
% (2) add standard normally distributed noise to it;
% (3) write out the data points.

function generate_data(k,b,npts)

% default parameters for the line
% k: slope; b: intercept
if nargin < 3
    k = 1;
    b = 1;
    npts = 50; % number of points
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create line

% x coord
x = linspace(1,50,npts);

% y coord
y = k * x + b;

% plot the line
figure(1);
plot(x,y,'k*-');
grid on;
axis equal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add noise

noise = randn(1,npts);
% Q: what if your data has a gaussian error, with mean of
% 0.1, and std of 0.1?
ynew = y + noise;

% add the new points to the figure
figure(1)
hold on;
plot(x,ynew,'bo');
xlabel('X'); ylabel('Y');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write out, simple matrix
dlmwrite('data.txt',[x(:) ynew(:) y(:)],'delimiter',' ');

% write out, with header
fid = fopen('data_h.txt','w+');
fprintf(fid,'%10s %10s %10s\n','xax','ydata','ytrue');
fprintf(fid,'%10.5f %10.5f %10.5f\n',[x(:) ynew(:) y(:)]');
fclose(fid);


