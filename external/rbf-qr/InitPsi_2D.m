function [Psi,C,Q,R,Z]=InitPsi_2D(ep,xk,tol)
% The purpose of this function is to compute \tilde{R} that is used
% for evaluating the basis functions \Psi used in the RBF-QR method.
%--- ep (scalar) : The shape parameter
%--- xk(1:N,1:2) : The center points in polar coordinates and
%                  scaled to the unit disc.
%figure(10),clf
  N = size(xk,1);  
  mp = eps; % machine precision
  tolD = 1e4*mp; % Regulates what is considered small in terms of pivoting out
  
%--- Find out the polynomial degree needed for including at least N columns
  jN=degree(N,2);

%--- Compute the first part of C, up to degree jN 
  [Psi,C,T] = AddCBlocks(ep,xk,jN);
  Ctot=C;
%--- QR-factorize the coefficient matrix and compute \tilde{R}
%--- incrementally block by block, until we find a block that is
%--- insignificant in absolute value compared to the previous ones.
%if (ep<=1)
%    tol = 2; % 100 times for pivoting
%  else
%    tol=5; % For large epsilon, sizes are more different and pivoting
%           % only needed for exactly singular columns
%  end  
  [Q,R]=IncQR(C,ep^2,Psi.j,tol);
  %  [R,Qfinal,columns,pow,extrapow]=testQR(C,ep,Psi)


  true=1;
  q = Q.q;
  posN=-1;
  jmax = jN;
  while true
    %--- Add one block to C
     jmax = jmax+1;
     [Psi,C,T] = AddCBlocks(ep,xk,jmax,Psi,T);
     Ctot = [Ctot C];

    %--- If we have less than N columns, then we cannot be done
    if (q > N)
      if (posN<0) % First time, find out location to compare with
	jN = Psi.j(R.order(N));
	pos = find(Psi.j==jN);
	posN = pos(1); % The location of d_j,0
      end		

      %--- Check the magnitude of the latest block. 
      jtest = Psi.j(end);
      pos = find(Psi.j==jtest);
      pos = pos(1);
      p1 = [1 posN]; % First and Nth, either one may be smallest
      p2 = [pos];
      relsc = EvalD_2D(ep,p1,p2,Psi.j,Psi.m,Psi.p,tolD);
      relsc = max(relsc);
      
      %--- Block small enough means we are done.
      if (relsc*exp(0.223*jtest+0.212-0.657*mod(jtest,2)) < 1e-16*mp)
	break
      end
    end  
    %--- Update the QR-factorization with the new block
    pos = find(Psi.j==Psi.j(end));
    [Q,R]=IncQR(C,ep^2,Psi.j(pos),tol,Q,R);
    q = Q.q;
  end
  M = length(R.order);
  % figure, H=pcolor(log10(abs(R.R))); axis ij
  %set(H,'EdgeColor','none'), colorbar
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Z=R.R(1:N,1:N)\Q.Q'*ones(N,1);
  p1 = R.order(1:N); p2=p1(end);
  D1 = EvalD_2D(ep,p1,p2,Psi.j,Psi.m,Psi.p,tolD);
  Z = diag(D1)*Z;
  Z = Z/max(abs(Z));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Rt = R.R(1:N,1:N)\R.R(1:N,N+1:M);
  p1 = R.order(1:N);    p2 = R.order((N+1):M);   
  if (M>N) % If there is a part to scale
    D = EvalD_2D(ep,p1,p2,Psi.j,Psi.m,Psi.p,tolD);
    %Dmax = max(max(abs(D)))
    Rt = D.*Rt;
  end
  %  D1 = 2*Psi.j(p1); % Extracting the powers
  D1 = EvalD_2D(ep,p1(1),p1,Psi.j,Psi.m,Psi.p,tolD);
  %Rtmax=max(max(abs(Rt)))
  %figure, H=pcolor(log10(abs(Rt))); axis ij
  %set(H,'EdgeColor','none'), colorbar
  %title(['Rt, ep=' num2str(ep)])
%  if (ep==0)
%    Rt=0*Rt;
%  end  
  Psi.ep = ep;
  Psi.xk = xk;
  Psi.Rt = Rt;
  Psi.columns = R.order;

function [Psi,C,T] = AddCBlocks(ep,xk,jmax,Psi,T)
% Compute new blocks of the coefficient matrix C up to degree jmax
%--- ep (scalar) : The shape parameter
%--- xk(1:N,1:2) : The center points as in InitPsi_2D
%--- jmax (scalar) : The degree to stop at
%--- Psi (struct) : Supplied the second or higher call
%--- T (struct) : Supplied the second or higher call

  if (nargin==3) % First call, initialize all
    j0 = 0;
    T.Pk = ones(size(xk,1),jmax+1); % Fill with ones up to the first end
    T.rscale = exp(-ep^2*xk(:,1).^2); % Row scaling of C = exp(-ep^2*r_k^2)
    Psi.j = zeros(0,1);
    Psi.m = zeros(0,1);
    Psi.p = zeros(0,1);
    Psi.cs = zeros(0,1);
  else % Psi has already been initialized, check previous end point
    j0 = Psi.j(end) + 1;
  end  
%--- Local indices to compute for in this call
  j = zeros(0,1); m=zeros(0,1); p=zeros(0,1); odd=mod(j0+1,2);
  for k=j0:jmax
    odd = abs(odd-1);
    p = [p; odd*ones(k+1,1)]; 	
    j = [j; k*ones(k+1,1)];
    q(1:2:k+1,1) = (0:(k-odd)/2)';
    q(2:2:k+1,1) = (abs(odd-1):(k-odd)/2)';
    m = [m;q(1:k+1)];
  end
  
%--- Fill in trigs and powers that will be reused later
  T.Hkc(:,max(1,j0):jmax) = cos(xk(:,2)*(max(1,j0):jmax));  
  T.Hks(:,max(1,j0):jmax) = sin(xk(:,2)*(max(1,j0):jmax));  
  for k=max(1,j0):jmax
    T.Pk(:,k+1) = xk(:,1).*T.Pk(:,k);
  end

  cs = zeros(size(j)); % Find which positions are sine and cosines
  pos = find(2*m+p>0); % cos=1, sin=-1
  cs(pos(1:2:end))=1;    cs(pos(2:2:end))=-1;

%--- Compute the new blocks of the coefficient matrix
  M = length(j); 
  cscale = 2*ones(1,M);           % Column scaling of C = b_{j,m}
  pos = find(2*m+p == 0);   cscale(pos) = 0.5*cscale(pos);
  pos = find(j-2*m == 0);   cscale(pos) = 0.5*cscale(pos);
  C = T.Pk(:,j+1); % The powers of r_k and then the trig part
  pos = find(cs ==  1);  C(:,pos) = C(:,pos).*T.Hkc(:,2*m(pos)+p(pos));
  pos = find(cs == -1);  C(:,pos) = C(:,pos).*T.Hks(:,2*m(pos)+p(pos));
  C = C.*(T.rscale*cscale);
  a = (j-2*m+p+1)/2; b=[j-2*m+1 (j+2*m+p+2)/2];
  z = ep.^4*xk(:,1).^2;
  for k=1:M
    C(:,k) = C(:,k).*hypergeom12(a(k),b(k,:),z);
  end  
  Psi.j = [Psi.j; j];
  Psi.m = [Psi.m; m];
  Psi.p = [Psi.p; p];
  Psi.cs = [Psi.cs; cs];

