clear, clc

% Matriz de exemplo
% A = [
%   1, 2, -1;
%   4, 3,  1;
%   2, 2,  3;
% ];

function [L,U] = decomposicao_LU(A)
  n = length(A);
  L = eye(n);
  U = A;

  for i=1:n-1
    for j=i+1:n
      alpha = U(j,i)/U(i,i);
      for k=i:n
        U(j,k) -= alpha*U(i,k);
      endfor

      L(j,i) = alpha;
      
      endfor
  endfor
  
endfunction

% Testando se isso deu certo
% [L,U] = decomposicao_LU(A);
% L, U