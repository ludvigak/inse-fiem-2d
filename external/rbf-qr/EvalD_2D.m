function D=EvalD_2D(ep,p1,p2,j,m,p,tol)
% The purpose of this function is to compute the scaling effect of 
% D_1^{-1} and D_2 applied to the correction matrix in the RBF-QR method.
%--- ep (scalar) : The shape parameter
%--- p1 (vector) : Indices for elements in D_1 
%--- p2 (vector) : Indices for elements in D_2
%--- j,m,p (vectors) : Identifiers for the expansion functions T_{j,m}  

%  tol=eps; % The smallest number we will agree to invert
  
  D = zeros(length(p1),length(p2));
  
  ep2 = 0.5*ep*ep;
  pmin = max(j(p1)); % Largest negative power that could occur
  if (ep2 < 1) % May need to limit the powers
    if (ep2 > tol)
      pmin = min(pmin, floor(log(tol)/log(ep2)));
    else    
      pmin=0;
    end
  end
%--- Precompute all positive and negative powers of ep that are present.
  pmax = max(j(p2));

  epp(0+pmin+1)=1;
  for pp=1:pmax
    epp(pp+pmin+1)=ep2*epp(pp+pmin);
  end
  ep2=1/ep2;
  for pp=1:pmin
    epp(-pp+pmin+1)=ep2*epp(-pp+pmin+2);
  end  

%--- Precompute powers of 2, both negative and postive
  mmax = max(m(p2));
  mmin = max(m(p1)); % Largest negative
  twop(0+mmin+1)=1;  % 2^0
  for pp=1:mmax
    twop(pp+mmin+1)=4*twop(pp+mmin);
  end
  for pp=1:mmin
    twop(-pp+mmin+1)=0.25*twop(-pp+mmin+2); 
  end

 %--- The column values stay the same for each row
  powj=j(p2);
  powm=m(p2);
  for k=1:length(p1) % For each row
    pow = powj-j(p1(k));
%--- Remove too negative powers
    pos = find(pow>=-pmin);
    D(k,pos) = epp(pow(pos)+pmin+1);
    pow = powm-m(p1(k));
    D(k,:) = D(k,:).*twop(pow+mmin+1);
  end
  
%--- This part does the ratios of factorials 
  f1 = (j(p2)+2*m(p2)+p(p2))/2;
  f2 = (j(p2)-2*m(p2)-p(p2))/2;
  f3 = (j(p1)+2*m(p1)+p(p1))/2;
  f4 = (j(p1)-2*m(p1)-p(p1))/2;
  fmax = max(max(f1),max(f3)); 
  fp = cumprod([1 1:fmax]); % Because 0!=1 
  for k=1:length(p1) % For each row
    numer = ones(1,size(D,2));
    denom = numer;
    pos = find(f3(k) < f1);
    denom(pos) = (1/fp(f3(k)+1))*fp(f1(pos)+1);
    pos = find(f4(k) < f2);
    denom(pos) = (1/fp(f4(k)+1))*denom(pos).*fp(f2(pos)+1);
    pos = find(f3(k) > f1);
    numer(pos) = fp(f3(k)+1)./fp(f1(pos)+1);
    pos = find(f4(k) > f2);
    numer(pos) = numer(pos).*(fp(f4(k)+1)./fp(f2(pos)+1));
    D(k,:) = D(k,:).*numer./denom;
  end

  


