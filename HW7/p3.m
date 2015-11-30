clear all; clc;
%% This is the solution script to problem 3 on homework set 3

%% Problem part (a)
% Define the input data and Tikhonov parameters
G1 = [1 1; 1 -3; 1 4; 1 5];
G2 = [1 -0.1; 1 0.3; 1 -0.4; 1 0.5];
G3 = [1 101; 1 97; 1 104; 1 105];
G4 = [1 3; 1 3; 1 3; 1 3];

alpha1 = 0.1;
alpha2 = 0.5;

% Calculate the regularized inverse matrices
reg1_1 = inv(G1'*G1 + alpha1^2*eye(2,2))*G1';
reg1_2 = inv(G2'*G2 + alpha1^2*eye(2,2))*G2';
reg1_3 = inv(G3'*G3 + alpha1^2*eye(2,2))*G3';
reg1_4 = inv(G4'*G4 + alpha1^2*eye(2,2))*G4';
reg2_1 = inv(G1'*G1 + alpha2^2*eye(2,2))*G1';
reg2_2 = inv(G2'*G2 + alpha2^2*eye(2,2))*G2';
reg2_3 = inv(G3'*G3 + alpha2^2*eye(2,2))*G3';
reg2_4 = inv(G4'*G4 + alpha2^2*eye(2,2))*G4';

%% Problem part (c)
%Define the data sets
d1 = [10 11 9 12]';
d2 = [10.1 11.4 8.7 9.8]';

%The 8 previous model parameter vectors are
m_d1_1  = [10.523 -0.013]'; 
m_d1_1t = [0.217 1.489]';
m_d2_1  = [10.465 -0.266]';
m_d2_1t = [0.180 1.234]';
m_d1_3  = [11.813 -0.013]';
m_d1_3t = [0.001 0.103]';
m_d2_3  = [37.046 -0.266]';
m_d2_3t = [9.616e-4 9.793e-2]';

%Calculate the regularized inverse for the two data sets and G1/G3
m_d1_1_r1 = reg1_1*d1;
m_d1_1_r2 = reg2_1*d1;
m_d2_1_r1 = reg1_1*d2;
m_d2_1_r2 = reg2_1*d2;
m_d1_3_r1 = reg1_3*d1;
m_d1_3_r2 = reg2_3*d1;
m_d2_3_r1 = reg1_3*d2;
m_d2_3_r2 = reg2_3*d2;

%Now set up the four different plots for the different combinations of 
%data sets and matrices

%For G1 using d1
x = linspace(-10,10,100);
y_d1_1_1 = m_d1_1(1) + m_d1_1(2)*x;
y_d1_1_1t = m_d1_1t(1) + m_d1_1t(2)*x;
y_d1_1_r1 = m_d1_1_r1(1) + m_d1_1_r1(2)*x;
y_d1_1_r2 = m_d1_1_r2(1) + m_d1_1_r2(2)*x;
x_d1 = G1(:,2);

plot(x_d1, d1, 'ko', ...
     x, y_d1_1_1, 'b-', ...
     x, y_d1_1_1t, 'r-', ...
     x, y_d1_1_r1, 'c-', ...
     x, y_d1_1_r2, 'm-', ...
     x_d1, d1, 'k.', ...
     [-1 1]*1000, [0 0], 'k--', ...
     [0 0], [-1 1]*1000, 'k--', 'MarkerSize', 10, 'LineWidth', 2);
legend('data', 'model_{generalized}', 'model_{truncated}', ...
      'model_{\alpha=0.1}', 'model_{\alpha=0.5}', 'Location', 'SouthEast');
title('matrix G_1 with data d_1');
xlabel('x');
ylabel('y');
box on;
axis([-10 10 -15 15])

print('-deps2c','-painters', 'p3b1');

%For G1 using d1
x = linspace(-10,10,100);
y_d2_1_1 = m_d2_1(1) + m_d2_1(2)*x;
y_d2_1_1t = m_d2_1t(1) + m_d2_1t(2)*x;
y_d2_1_r1 = m_d2_1_r1(1) + m_d2_1_r1(2)*x;
y_d2_1_r2 = m_d2_1_r2(1) + m_d2_1_r2(2)*x;
x_d2 = G1(:,2);

plot(x_d2, d2, 'ko', ...
     x, y_d2_1_1, 'b-', ...
     x, y_d2_1_1t, 'r-', ...
     x, y_d2_1_r1, 'c-', ...
     x, y_d2_1_r2, 'm-', ...
     x_d2, d2, 'k.', ...
     [-1 1]*1000, [0 0], 'k--', ...
     [0 0], [-1 1]*1000, 'k--', 'MarkerSize', 10, 'LineWidth', 2);
legend('data', 'model_{generalized}', 'model_{truncated}', ...
      'model_{\alpha=0.1}', 'model_{\alpha=0.5}', 'Location', 'SouthEast');
title('matrix G_1 with data d_2');
xlabel('x');
ylabel('y');
box on;
axis([-10 10 -15 15])

print('-deps2c','-painters', 'p3b2');

%For G3 using d1
x = linspace(-10,120,100);
y_d1_3_1 = m_d1_3(1) + m_d1_3(2)*x;
y_d1_3_1t = m_d1_3t(1) + m_d1_3t(2)*x;
y_d1_3_r1 = m_d1_3_r1(1) + m_d1_3_r1(2)*x;
y_d1_3_r2 = m_d1_3_r2(1) + m_d1_3_r2(2)*x;
x_d1 = G3(:,2);

plot(x_d1, d1, 'ko', ...
     x, y_d1_3_1, 'b-', ...
     x, y_d1_3_1t, 'r-', ...
     x, y_d1_3_r1, 'c-', ...
     x, y_d1_3_r2, 'm-', ...
     x_d1, d1, 'k.', ...
     [-1 1]*1000, [0 0], 'k--', ...
     [0 0], [-1 1]*1000, 'k--', 'MarkerSize', 10, 'LineWidth', 2);
legend('data', 'model_{generalized}', 'model_{truncated}', ...
      'model_{\alpha=0.1}', 'model_{\alpha=0.5}', 'Location', 'SouthEast');
title('matrix G_3 with data d_1');
xlabel('x');
ylabel('y');
box on;
axis([-10 120 -10 20])

print('-deps2c','-painters', 'p3b3');

%For G3 using d1
x = linspace(-10,120,100);
y_d2_3_1 = m_d2_3(1) + m_d2_3(2)*x;
y_d2_3_1t = m_d2_3t(1) + m_d2_3t(2)*x;
y_d2_3_r1 = m_d2_3_r1(1) + m_d2_3_r1(2)*x;
y_d2_3_r2 = m_d2_3_r2(1) + m_d2_3_r2(2)*x;
x_d2 = G3(:,2);

plot(x_d2, d2, 'ko', ...
     x, y_d2_3_1, 'b-', ...
     x, y_d2_3_1t, 'r-', ...
     x, y_d2_3_r1, 'c-', ...
     x, y_d2_3_r2, 'm-', ...
     x_d2, d2, 'k.', ...
     [-1 1]*1000, [0 0], 'k--', ...
     [0 0], [-1 1]*1000, 'k--', 'MarkerSize', 10, 'LineWidth', 2);
legend('data', 'model_{generalized}', 'model_{truncated}', ...
      'model_{\alpha=0.1}', 'model_{\alpha=0.5}', 'Location', 'SouthEast');
title('matrix G_3 with data d_2');
xlabel('x');
ylabel('y');
box on;
axis([-10 120 -10 20])

print('-deps2c','-painters', 'p3b4');
