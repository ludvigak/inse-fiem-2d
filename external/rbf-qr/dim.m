function N=dim(K,n)    
%
% Find the dimension of the polynomial space of degree K in n dimensions
%  
  N = prod([(K+1):(K+n)]./[1:n]);
