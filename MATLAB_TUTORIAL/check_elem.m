% check the (i,j) element 
function flag = check_elem(i,j,val)
if i == j
    flag = val == 1; 
else
    flag = val == 0;
end
end