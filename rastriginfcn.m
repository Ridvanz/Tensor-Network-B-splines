% Computes the value of Rastrigin benchmark function.
% SCORES = RASTRIGINFCN(X) computes the value of the Rastrigin function at 
% point X. RASTRIGINFCN accepts a matrix of size M-by-N and returns a vetor 
% SCORES of size M-by-1 in which each row contains the function value for
% the corresponding row of X.
% For more information please visit: 
% https://en.wikipedia.org/wiki/Rastrigin_function
% 
% Author: Mazhar Ansari Ardeh
% Please forward any comments or bug reports to mazhar.ansari.ardeh at
% Google's e-mail service or feel free to kindly modify the repository.
function f = rastriginfcn(x)
    x = (x-0.5)*10;
    n = size(x, 2);
    A = 10;
    f = ((A * n) + (sum(x .^2 - A * cos(2 * pi * x), 2)))/100;
end