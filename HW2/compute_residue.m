function [ r ] = compute_residue( x,y,M,d )

xs = M(1);
ys = M(2);
zs = M(3);
p = M(4);
r =d- p*zs./((x - xs).^2 + (y - ys).^2 + zs^2).^(3/2);
end

