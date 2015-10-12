function [root,flag] = find_root(x0,epsilon,itermax)
root = x0;
flag = false;   % whether you find a root
for i = 1:itermax
    f_curr = f_func(x0);
    f_grad_curr = f_grad_func(x0);
    
    if abs(f_grad_curr) < epsilon   % denominator is too small, stop
       break; 
    end
    
    x1 = x0 - f_curr/f_grad_curr;   
    
    if abs( x1-x0 ) < epsilon
        flag = true;     % find the root
        root = x0;
        break;
    end
    
    x0 = x1;
end
end

% value of the function at x
function val = f_func(x)
val = x^2/10 + 6 * sin(x) + 2;
end

% value of the gradient at x
function val = f_grad_func(x)
val = x/5 + 6 * cos(x);
end