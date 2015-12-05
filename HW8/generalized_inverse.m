function [model,Resolution] = generalized_inverse(data,GG,p)

if nargin <= 2  % if p is not set, do not truncate
    p = min(size(GG));
end

[UU,DD,VV] = svd(GG);

I_DD = zeros(size(DD'));

for tmp = 1:1:p 
    
I_DD(tmp,tmp) = 1/DD(tmp,tmp); 

end


model = VV*I_DD*UU'*data;

Resolution = VV*I_DD*UU'*UU*DD*VV';

end


