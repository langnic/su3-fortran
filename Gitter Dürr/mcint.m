function[] = mcint(n)

  if nargin < 1, n = 1000; end
  
  f = @(t) t.^(-1/3) + t./10;
  x = rand(n,1);
  f1 = mean(f(x));
  f2 = mean(f(x).^2);
  sigma = sqrt(f2 - f1^2)/sqrt(n);
  
  g = @(t) 3/2 + 3/20 * t.^2;
  x = rand(n,1);
  g1 = mean(g(x));
  g2 = mean(g(x).^2);
  sigma = sqrt(g2 - g1^2)/sqrt(n-1)
  
end


