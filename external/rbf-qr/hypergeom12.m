function y=hypergeom12(a,b,x)
%
% In our case, I will only implement the 1F2 case, which we have.
%
  if (length(a)~=1 | length(b)~=2)
    error('Wrong number of arguments in hypergeom for 1F2')
  end
  %
  % The first coefficient is 1 always
  % 
  alpha=1;
  [ma,pos]=max(abs(x)); % Could possibly be complex
  v = ones(size(x));
  y = v;
  n=0;
  test=1;
  mp = eps; % machine precision
  while (abs(test) > mp)
    %
    % The hypergeometric coefficient is computed from the arguments.
    % http://en.wikipedia.org/wiki/Hypergeometric_function#The_series_pFq
    %
    alpha = alpha*(a+n)/(b(1)+n)/(b(2)+n);
    n = n+1;
    %
    % Power and factorial can also be done recursively.
    %
    v = (1/n)*v.*x;
    y = y + alpha*v;
    test = alpha*v(pos);
  end
  
    