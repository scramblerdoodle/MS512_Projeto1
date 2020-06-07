clc,clear
for i = 0:20
   x(i+1)=i*0.05;
   b(i+1,1)=exp(sin(6*i*0.05));
endfor
A=vander(x,11);
%Método 1 Cholesky
R1=chol(A'*A);
opts.LT=true;
y=linsolve(R1',A'*b,opts);
opts2.UT=true;
c1=linsolve(R1,y,opts2);

%Método 2 QR condensada
V=A'*A;
for j=1:11
  for i=1:j-1
    R2(i,j)=Q1(:,i)'*V(:,j);
    V(:,j)=V(:,j)-R2(i,j)*Q1(:,i);
  endfor
  R2(j,j)=norm(V(:,j));
  if(R2(j,j)==0)
    disp("Deu erro");
    break;
  else
    Q1(:,j)=V(:,j)/R2(j,j);
  endif
endfor
c2=linsolve(R2, Q1'*A'*b,opts2);

%Método 3 QR completa
[Q,R3] = qr(A'*A);
c3=linsolve(R3, Q'*A'*b,opts2);

%Comparação
a1=0.5*norm(A*c1-b)^2
a2=0.5*norm(A*c2-b)^2
a3=0.5*norm(A*c3-b)^2