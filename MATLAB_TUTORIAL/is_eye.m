function is_eye(a_mat)

% Default a_mat, if no input
if nargin < 1
    a_mat = magic(3);
    %a_mat = zeros(2,3);
end

% Check size first
% if/else
if size(a_mat,1) ~= size(a_mat,2)
    disp('not an identity matrix!');
    return;
else
    disp('size is OK!');
    N = size(a_mat,1);
end

% Check element
% for loop
for i = 1:N
    for j = 1:N
        flag = check_elem(i,j,a_mat(i,j));
        if ~flag
            disp('NOT an identity matrix!');
            return
        end
    end
end

disp('IS an identity matrix!');  
disp(a_mat);

% Just for illustration of
% while loop
i = 1;
while i <= min(100,N)
    disp(a_mat(i,i));
    i = i+1;
%     if i > 100      % will only display the first 100 diag elements
%         break
%     end
end
end



    