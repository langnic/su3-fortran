function [a,b] = miesman(A,b)
a = det(A);
b = A\b;
end
