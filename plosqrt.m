function S = plosqrt(A)
%POLSQRT Matrix square root
%   S = plosqrt(A) Computes the square root S = A^(1/2) of the symmetric
%   positive definite matrix A by doing a Cholesky decomposition and
%   calling poldec.

[R, p] = chol(A);
if p~= 0
    disp('sorry, I only work for a sym.pos.def matrix.');
else
    [~, S] = poldec(R);
    clc; % clear the displays that result from calling poldec()
    disp('>> plosqrt()');
end
end

