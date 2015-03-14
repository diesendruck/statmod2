[t1sq,b] = meshgrid(-1:.1:1);
d <- sqrt(sum((i2 - i1)^2))
val <- t1.sq*(1 + sqrt(5)*d/b + 5/3*(d^2)/(b^2)) * exp(-sqrt(5)*d/b) + (t2.sq)*Delta(i2, i1)

n1 = (1:1:10);
n2 = (0:.1:1);
z = zeros(length(n1),length(n2));
for i = 1:n1
  for j = 1:n2
    d = sqrt(sum((i-j)^2));
    z(i,j) = d;
  end
end


surf(t1sq,b,z);
