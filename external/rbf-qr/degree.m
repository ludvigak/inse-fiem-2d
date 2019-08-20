function K=degree(N,n)
%
% Find the polynomial degree that N basis functions correspond to in n
% dimensions
%  
  for k=0:N-1 % K(N) cannot be larger than N-1 (1D-case)
    if (dim(k,n) >= N)
      K=k;
      break
    end
  end
  
