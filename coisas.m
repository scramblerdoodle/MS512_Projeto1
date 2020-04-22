1;

function [ENV, ENVcol, ENVlin] = env(A)
  %ENVELOPEcoluna ENVELOPEcoluna ENVELOPEcoluna
  n = length(A);
  aux = 0;
  count = 0;
  first = 0;
  charge = 1;
  auxLin = 0;
  ENV = ENVlin = 0;
  ENVcol = [1 1];

  %ENV = zeros(1,floor((n*n)/2-n));
  ENVcol = zeros(1,n+1);
  %ENVlin = zeros(1,floor((n*n/2)-n));


  for j = 2:n
    for i = 1:j-1
      if (A(i,j) != 0)
        auxLin = i;
        break;
      endif
    endfor
    if (auxLin != 0)
      for k = auxLin:j-1
        count++;
        if(k == auxLin)
        first = count;
        endif
        ENV(count) = A(k,j); 
        ENVlin(count) = k;
      endfor
      
      ENVcol(j) = first;
      for k = j+1:n+1
        ENVcol(k) = j-auxLin+ENVcol(j);
      endfor
      while(charge>0)
        ENVcol(j-charge) = first;
        charge--;
      endwhile
    else
      charge++;
    endif
    auxLin = 0;
  endfor
  %ENVELOPEcoluna ENVELOPEcoluna ENVELOPEcoluna
endfunction
  
  
function [ENVi, ENVlini, ENVcoli, ENVs, ENVcols, ENVlins, DIAG] = envelope(A)
  [ENVi, ENVlini, ENVcoli] = env(A');
  [ENVs, ENVcols, ENVlins] = env(A);
  for i = 1:length(A)
    DIAG(i) = A(i,i);
  endfor
endfunction
  

% TODO: essa função não tá sendo usada
function [x] = sistema_linearCol(U,b)
  %SISTEMA LINEAR SISTEMA LINEAR
  [ENV, ENVcol, ENVlin] = env(U);
  n = length(U);
  x = zeros(n,1);
  for i = n:-1:1
    x(i) = b(i)/U(i,i);
    if ENVcol(i)!=ENVcol(i+1)
      for j = ENVcol(i):ENVcol(i+1)-1
        b(ENVlin(j)) -=  x(i)*ENV(j);
      endfor
    endif
  endfor
  %SISTEMA LINEAR SISTEMA LINEAR
endfunction


function [x] = sistema_linearLin(ENV, ENVlin, ENVcol, DIAG, b)
  %SISTEMA LINEAR SISTEMA LINEAR
  n = length(DIAG);
  x = zeros(n,1);
  for i = 1:n 
    s = b(i);
    if ENVlin(i)!=ENVlin(i+1)
      for j = ENVlin(i):ENVlin(i+1)-1 
        s -=  x(ENVcol(j))*ENV(j);
      endfor
    endif
    x(i) = s/DIAG(i);
  endfor
  %SISTEMA LINEAR SISTEMA LINEAR
endfunction


function [b,c] = theFinder(ENVi, ENVlini, ENVcoli, ENVs, ENVcols, ENVlins, i);
  k = i;
  b = zeros(i);
  c = zeros(i);
  for j = ENVcols(i+2)-1:-1:ENVcols(i+1)
    b(k) = ENVs(j);
    k = k-1;  
  endfor
  k = i;
  for j = ENVlini(i+2)-1:-1:ENVlini(i+1)
    c(k) = ENVi(j);
    k = k-1; 
  endfor  
endfunction


function [ENV,ENVcol,ENVlin] = newColumn(ENV, ENVcol, ENVlin, b)
  i = 1;
  while(i <= length(b) && b(i) == 0  ) 
    i = i+1;
  endwhile
  aux = length(ENVcol);
  for j = 0:length(b)-i
    ENV(ENVcol(aux)+j) = b(i+j);
    ENVlin(ENVcol(aux)+j) = i+j;
  endfor
  ENVcol(aux+1) = ENVcol(aux)+length(b)-i+1;
endfunction 


function [ENVL,ENVlinL,ENVcolL,DIAGL,ENVU,ENVcolU,ENVlinU,DIAGU] = LUENVELOPE(ENViA, ENVliniA, ENVcoliA, ENVsA, ENVcolsA, ENVlinsA, DIAGA)
  ENVlinL = [1 1];
  n = length(DIAGA);
  DIAGL(1) = 1;
  DIAGU(1) = DIAGA(1);
  ENVcolU = [1 1];
  ENVL = 0;
  ENVcolL = 0;
  ENVU = 0; 
  ENVlinU = 0;
  ENVcolU = [1 1];
  
  for i = 1:n-1
    [b,c] = theFinder(ENViA, ENVliniA, ENVcoliA, ENVsA, ENVcolsA, ENVlinsA, i);
    u = sistema_linearLin(ENVL, ENVlinL, ENVcolL,DIAGL,b);
    l = sistema_linearLin(ENVU, ENVcolU, ENVlinU,DIAGU,c);
    [ENVL,ENVlinL,ENVcolL] = newColumn(ENVL,ENVlinL,ENVcolL,l);
    [ENVU,ENVcolU,ENVlinU] = newColumn(ENVU,ENVcolU,ENVlinU,u);
    DIAGL(i+1) = 1;
    s = 0;
    for k = 1:length(l)
      s = s + l(k)*u(k);
    endfor
    DIAGU(i+1) = DIAGA(i+1)-s;
  endfor

endfunction


% essa função é pra quê? hsdfhksd
function M = sparsity(U)
  n = length(U);
  for i = 1:n
    for j = 1:n
      if (U(i,j) == 0)
        M(i,j) = 0;
      else
        M(i,j) = 1;
      endif
    endfor
  endfor
endfunction


% LU por pivoteamento (?)
function [L,U,P] = lupp(A)
  n = length(A);
  p = (1:n)';
  for k = 1:n-1
    [temp,pos] = max(abs(A(k:n,k)));
    row2swap = k-1+pos;
    A([row2swap, k],:) = A([k, row2swap],:);
    p([row2swap, k]) = p([k, row2swap]);
    J = k+1:n;
    A(J,k) = A(J,k)/A(k,k);
    A(J,J) = A(J,J) - A(J,k)*A(k,J);
  endfor
  
  for i = 1:n
    P(i,p(i)) = 1;
  endfor
  
  U = triu(A);
  L = eye(length(A))+tril(A,-1);
endfunction


% LU sem pivoteamento
function [L,U] = lusp(M)
  n = length(M);
  U = M;
  L = eye(n);
  for j = 1:n-1
    for i = j+1:n
      L(i,j) = U(i,j)/U(j,j);
      U(i,j:n) = U(i,j:n) - L(i,j)*U(j,j:n);
    endfor
  endfor
endfunction
