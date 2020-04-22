1;

%{
  Funcao que monta o envelope da parte superior de uma matriz A por colunas
  
  INPUT: Matriz A
  OUTPUT: ENV, ENVcol e ENVlin de A
  
  (note que trocando os papeis de ENVcol e ENVlin, e aplicando em A transposto,
  temos o envelope da parte inferior de uma matriz A por linhas)
%}
function [ENV, ENVcol, ENVlin] = env(A)
  n = length(A);
  auxCol = 0; % Variavel que guarda o primeiro numero nao-nulo de uma coluna
  count = 0;  % Contador para os indices usados em ENV e ENVlin
  first = 0;  % Indice do primeiro numero nao-nulo em uma coluna
  charge = 1; % "Carga" que acumula cada vez que temos uma coluna nula, pois os indices repetem em ENVcol nesse caso
  
  ENV = ENVlin = 0;
  ENVcol = zeros(1,n+1);


  for j = 2:n

    % Encontra o primeiro numero nao-nulo da parte superior da coluna j
    for i = 1:j-1
      if (A(i,j) != 0)
        auxCol = i;
        break;
      endif
    endfor

    % Se for uma linha nao-nula, salva todos os valores (a partir do primeiro numero nao-nulo) no envelope
    if (auxCol != 0)
      for k = auxCol:j-1
        count++;
        if(k == auxCol)
          first = count;
        endif
        ENV(count) = A(k,j); 
        ENVlin(count) = k;
      endfor
      
      % Monta o ENVcol
      ENVcol(j) = first;
      for k = j+1:n+1
        ENVcol(k) = j-auxCol+ENVcol(j);
      endfor

      while(charge>0)
        ENVcol(j-charge) = first;
        charge--;
      endwhile

    else
      charge++;

    endif
    
    % Reseta auxCol para a proxima passagem
    auxCol = 0;
  endfor
endfunction


%{ 
  Funcao que monta os envelopes inferiores, envelopes superiores e a diagonal de A

  INPUT: Matriz A
  OUTPUT: ENV_inferior, ENVlin_inferior, ENVcol_inferior,
          ENV_superior, ENVcol_superior, ENVlin_superior,
          DIAG
%}
function [ENVi, ENVlini, ENVcoli, ENVs, ENVcols, ENVlins, DIAG] = envelope(A)
  [ENVi, ENVlini, ENVcoli] = env(A');
  [ENVs, ENVcols, ENVlins] = env(A);
  for i = 1:length(A)
    DIAG(i) = A(i,i);
  endfor
endfunction
  

% Solucao de sistema linear por colunas a partir de envelopes
function [x] = sistema_linearCol(ENV, ENVcol, ENVlin, DIAG, b)
  n = length(DIAG);
  x = zeros(n,1);
  for i = n:-1:1
    x(i) = b(i)/DIAG(i);
    if ENVcol(i) != ENVcol(i+1)
      for j = ENVcol(i):ENVcol(i+1)-1
        b(ENVlin(j)) -=  x(i)*ENV(j);
      endfor
    endif
  endfor
endfunction


% Solucao de sistema linear por linhas a partir de envelopes
function [x] = sistema_linearLin(ENV, ENVlin, ENVcol, DIAG, b)
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
endfunction


%{
  Encontra os 'subvetores' b e c na decomposicao explicada no exercicio 5

  INPUT: Envelopes inferiores e superiores de uma matriz, e o indice i que estamos analisando
  OUTPUT: Vetores b e c, usados em LUENVELOPE
%}
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


% Adiciona uma nova coluna b ao envelope
function [ENV,ENVcol,ENVlin] = newColumn(ENV, ENVcol, ENVlin, b)
  i = 1;
  while(i <= length(b) && b(i) == 0) 
    i = i+1;
  endwhile
  aux = length(ENVcol);
  for j = 0:length(b)-i
    ENV(ENVcol(aux)+j) = b(i+j);
    ENVlin(ENVcol(aux)+j) = i+j;
  endfor
  ENVcol(aux+1) = ENVcol(aux)+length(b)-i+1;
endfunction 

%{
  Algoritmo que realiza todo o envelope da decomposicao LU a partir dos envelopes de uma matriz A

  INPUT:  ENV_inferior_A, ENVlin_inferior_A, ENVcol_inferior_A,
          ENV_superior_A, ENVcol_superior_A, ENVlin_superior_A,
          DIAG_A

  OUTPUT: todos os envelopes de L e de U
    i.e.  ENV_L, ENVlin_L, ENVcol_L, DIAG_L,
          ENV_U, ENVcol_U, ENVlin_U, DIAG_U 

%}
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


% essa funcao é pra quê? hsdfhksd
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
