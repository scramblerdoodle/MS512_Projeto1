A = [1 0 0 0 0 0 0;
     0 1 0 0 1 0 0;
     0 0 1 0 0 0 0;
     0 1 0 1 0 1 0;
     0 0 0 0 1 0 0;
     0 1 0 0 0 1 1;
     0 0 0 0 0 0 1;]
     
function [ENV, ENVcol, ENVlin] = envelope(U)
  %ENVELOPE ENVELOPE ENVELOPE
  n = length(U);
  aux = 0;
  count = 0;
  first = 0;
  charge = 1;
  auxLin = 0;

  ENV = zeros(1,floor((n*n)/2-n));
  ENVcol = zeros(1,n+1);
  ENVlin = zeros(1,floor((n*n/2)-n));


  for j = 2:n
    for i = 1:j-1
      if (U(i,j) != 0)
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
        ENV(count) = U(k,j); 
        ENVlin(count) = k;
      endfor
      
      ENVcol(j) = first;
      while(charge>0)
        ENVcol(j-charge) = first;
        charge--;
      endwhile
    else
      charge++;
    endif
    if(j!=n)
      auxLin = 0;
    else
      ENVcol(j+1) = ENVcol(j) + (j-auxLin);
    endif
  endfor
  %ENVELOPE ENVELOPE ENVELOPE
endfunction