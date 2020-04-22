n = 13;
A=zeros(n);
vetor_perm = [2 3 1 4 5 7 6 9 10 11 8 12 13];

function A = forcas(A)
  alpha = 1/sqrt(2);
  # i = número da equação
  # j = variável associada
  
  # f2 - f6 = 0
  A(1,2) = 1;
  A(1,6) = -1;
  
  # f3 = 10
  A(2,3) = 1;
  
  # alpha*f1 - f4 - alpha*f5 = 0
  A(3,1) = alpha;
  A(3,4) = -1;
  A(3,5) = -alpha;
  
  # alpha*f1 + f3 + alpha*f5 = 0
  A(4,1) = alpha;
  A(4,3) = 1;
  A(4,5) = alpha;
  
  # f4 - f8 = 0
  A(5,4) = 1;
  A(5,8) = -1;
  
  # f7 = 0
  A(6,7) = 1;
  
  # alpha*f5 + f6 - alpha*f9 - f10 = 0
  A(7,5) = alpha;
  A(7,6) = 1;
  A(7,9) = -1;
  A(7,10) = -1;
  
  # alpha*f5 + f7 + alpha*f9 = 15
  A(8,5) = alpha;
  A(8,7) = 1;
  A(8,9) = alpha;
  
  # f10 - f13 = 0
  A(9,10) = 1;
  A(9,13) = -1;
  
  # f11 = 20
  A(10,11) = 1;
  
  # f8 + alpha*f9 - alpha*f12 = 0
  A(11,8) = 1;
  A(11,9) = alpha;
  A(11,12) = -1;
  
  # alpha*f9 + f11 + alpha*f12 = 0
  A(12,9) = alpha;
  A(12,11) = 1;
  A(12,12) = 1;
  
  # f13 - alpha*f12 = 0
  A(13,13) = 1;
  A(13,12) = -alpha;
endfunction

function M = constroi_por_envelope(DIAG, ENV, ENVcol, ENVlin)
  M = zeros(length(DIAG));
   for j = 1:length(ENVcol)-1
      for k = ENVcol(j): ENVcol(j+1)-1
    M(ENVlin(k), j) = ENV(k);  
    endfor
  endfor
  
  for i = 1:length(DIAG)
    M(i,i) = DIAG(i);
  endfor
endfunction

A = forcas(A);
# spy(A);
P = eye(n) (:, vetor_perm);
# spy(P*A);

proj1;
[ENVi, ENVlini, ENVcoli, ENVs, ENVcols, ENVlins, DIAG] = envelope(P*A)
[ENVL,ENVlinL,ENVcolL,DIAGL,ENVU,ENVcolU,ENVlinU,DIAGU] = LUENVELOPE(ENVi, ENVlini, ENVcoli, ENVs, ENVcols, ENVlins,DIAG)


U = constroi_por_envelope(DIAGU, ENVU, ENVcolU, ENVlinU)
# spy(U);

L = constroi_por_envelope(DIAGL, ENVL, ENVlinL, ENVcolL)'
# spy(L);
