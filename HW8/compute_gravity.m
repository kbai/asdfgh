function [ gz ] = compute_gravity( m,D,x )
  %number of measurements
  gz = zeros(size(x));    
  G=6.67e-11;  %gravity constant
  L = length(m); %length of mass anormaly  equals 101 in our case
for ii = 1:1:L

    gz = gz + G*D*m(ii)./(D^2+(x-0.1*(ii-1)).^2).^(1.5);
end

end

