kernel = 1
switch kernel
  case 1; k = @(x,y) exp(-0.01*sqrt((x-y)'*(x-y))); %Squared Exponential b=sqrt(50)
  case 2; k = @(x,y) exp(-0.1*sqrt((x-y)'*(x-y))); %Squared Exponential b=sqrt(5)
  case 3; k = @(x,y) exp(-1*sqrt((x-y)'*(x-y))); %Squared Exponential b=1/sqrt(2)
  case 4; k = @(x,y) exp(-100*sqrt((x-y)'*(x-y))); %Squared Exponential b=1/sqrt(200)
  case 5; k = @(x,y) exp(-100*(x-y)'*(x-y)); %Laplace
  case 6; k = @(x,y) (1 + sqrt(5)*sqrt((x-y)'*(x-y))/.05 + 5/3*(sqrt((x-y)'*(x-y))^2)/(.05^1)) * exp(-sqrt(5)*sqrt((x-y)'*(x-y))/.05); %Matern b=0.05
  case 7; k = @(x,y) (1 + sqrt(5)*sqrt((x-y)'*(x-y))/.5 + 5/3*(sqrt((x-y)'*(x-y))^2)/(.5^1)) * exp(-sqrt(5)*sqrt((x-y)'*(x-y))/.5); %Matern b=0.5
  case 8; k = @(x,y) (1 + sqrt(5)*sqrt((x-y)'*(x-y))/3 + 5/3*(sqrt((x-y)'*(x-y))^2)/(3^2)) * exp(-sqrt(5)*sqrt((x-y)'*(x-y))/3); %Matern b=3
end

points = (0:0.05:1)';
[U, V] = meshgrid(points, points);
x = [U(:) V(:)]';
n = size(x,2);

C = zeros(n,n);
for i = 1:n
  for j = 1:n
    C(i,j) = k(x(:,i),x(:,j));
  end
end

u = randn(n,1);
[A,S,B] = svd(C);
z = A*sqrt(S)*u;

figure(0); clf;
Z = reshape(z, sqrt(n),sqrt(n));
surf(U,V,Z);