function [T,F]=RBF_QR_precomp_2D(j,m,p,ep,xe,deg,T)
%--- Extract descriptors for the Psi-functions
  tol = 10*eps; % 10 times the machine precision
  jmax = j(end);    pmax = p(end);
  Ne = size(xe,1);
  M = length(j);
  
  if (nargin <6 || isempty(T)) % First time
    %--- First compute the basic functions that T are built from
    T.F=[]; T.Q=[]; % Used for testing later
    T.deg = 0;
    T.Te = cos(acos(xe(:,1))*(0:jmax));
    
    T.Hec = cos(xe(:,2)*(1:jmax));
    T.Hes = sin(xe(:,2)*(1:jmax));
  
    T.xe = xe;
    T.re2 = xe(:,1).^2;  % Only even powers are needed for evaluation points
    T.Pe = ones(Ne,(jmax-pmax)/2+1);
    for pp=1:(jmax-pmax)/2
      T.Pe(:,pp+1) = T.re2.*T.Pe(:,pp);
    end  
    T.rsc = exp(-ep^2*T.re2);
    T.r = xe(:,1);
  end
  
  if (deg>=1 & T.deg<1)
    T.deg = 1;
    T.fac  = 1-T.re2;     
    pos = find(abs(T.r-1)<=tol); T.fac(pos)=1; T.fac=1./T.fac;

    T.dTe = -((T.r.*T.fac)*(0:jmax)).*T.Te;
    T.dTe(:,2:end) = T.dTe(:,2:end) + (T.fac*(1:jmax)).*T.Te(:,1:end-1);
    T.dTe(pos,:) = ones(length(pos),1)*(0:jmax).^2;
  end

  if (deg>=2 & T.deg<2)
    T.deg = 2;
    jj = (0:jmax);
    pos = find(abs(T.r-1)<=tol); 

    T.d2Te = -T.fac*jj.^2.*T.Te+(xe(:,1).*T.fac)*ones(1,length(jj)).*T.dTe;
    T.d2Te(pos,:) = ones(length(pos),1)*(jj.^2.*(jj.^2-1)/3);
  end
      
  if (deg==1 & isempty(T.F)) % Needed now, but not computed yet
    pr0 = find(T.r<=tol);
    pm0 = find(m==0); % For s=0, the coeff is zero
    pm = find(m>0);

    %--- (2m+p)T_{j,m}/r=(2m+p)exp(-ep^2r^2)r^{2m-1}T_{j-2m}(r)    
    T.G = zeros(Ne,M);
    
    T.G(:,pm) = (T.rsc.*T.r)*ones(1,length(pm)).* ...
	      T.Pe(:,m(pm)+1-1).* ( ... % Corresponds to 2*m-2
		  (-2*ep^2*T.re2*ones(1,length(pm))+2*ones(Ne,1)* ...
		   m(pm)').*T.Te(:,j(pm)-2*m(pm)+1)+T.r*ones(1, ...
		length(pm)).*T.dTe(:,j(pm)-2*m(pm)+1));
 
    T.G(:,pm0) = T.rsc*ones(1,length(pm0)).* ...
	( -2*ep^2*T.r*ones(1,length(pm0)).*T.Te(:,j(pm0)+1)+ ...
	  T.dTe(:,j(pm0)+1));

    T.F = zeros(Ne,M);
    T.F(:,pm) = (T.rsc.*T.r)*(2*m(pm)+p(pm))'.* ...
	      T.Pe(:,m(pm)+1-1).*T.Te(:,j(pm)-2*m(pm)+1);

    T.F(:,pm0) = (T.rsc./T.r)*(2*m(pm0)+p(pm0))'.*... % Either one or zero
	T.Te(:,j(pm0)+1);

    T.F(pr0,pm0) = ones(length(pr0),1)*(j(pm0).*cos((j(pm0)-1)/2*pi))';
  end
  
  if (deg==2 & isempty(T.Q))    
    pr0 = find(T.r<=tol);
    pm0 = find(m==0); % For s=0, the coeff is zero
    pm = find(m>0);

    T.Q = zeros(Ne,M);
    cT = 2*ep^2*(2*ep^2*T.re2*ones(1,M)-4*ones(Ne,1)*m' - 1 );
    cT(:,pm) = cT(:,pm).*(T.re2*ones(1,length(pm))) + ...
	ones(Ne,1)*(2*m(pm).*(2*m(pm)-1))';
    %
    cdT = -4*ep^2*T.r*ones(1,M);
    cdT(:,pm) = cdT(:,pm).*(T.re2*ones(1,length(pm)))+4*T.r*m(pm)';
    
    cd2T = ones(Ne,M);
    cd2T(:,pm) = cd2T(:,pm).*(T.re2*ones(1,length(pm)));
    
    T.Q = cT.*T.Te(:,j-2*m+1) + cdT.*T.dTe(:,j-2*m+1) + cd2T.*T.d2Te(:,j-2*m+1);
    T.Q = (T.rsc*ones(1,M)).*T.Q;
    T.Q(:,pm) = T.Pe(:,m(pm)-1 +1).*T.Q(:,pm);
      
    %--- The rest of the terms have special cases also for m=1
    rinv=T.r; rinv(pr0)=1; rinv=1./rinv; % Max value is 14.
    T.R = zeros(Ne,M);
      
    cT = ones(Ne,1)*(1-2*m)';
    cT(:,pm) = cT(:,pm) + 2*ep^2*(T.re2*ones(1,length(pm)));
    cT(:,pm0) = cT(:,pm0).*(rinv.^2*ones(1,length(pm0))) + 2*ep^2;
      
    cdT = zeros(Ne,M);
    cdT(:,pm) = -T.r*ones(1,length(pm));
    cdT(:,pm0) = -rinv*ones(1,length(pm0));

    T.R = cT.*T.Te(:,j-2*m+1) + cdT.*T.dTe(:,j-2*m+1);
    T.R(pr0,pm0) = 0; % Limit for m=r=0 is zero in all cases
    T.R = (T.rsc*(2*m+p)').*T.R; % This makes the m=p=0 case = 0
    T.R(:,pm) = T.R(:,pm).*T.Pe(:,m(pm)-1 +1);
          
    T.S = zeros(Ne,M);
      
    cT = ones(Ne,1)*(2*m-(2*m+p).^2)';
    cT(:,pm) = cT(:,pm)-2*ep^2*T.re2*ones(1,length(pm));
    cT(:,pm0) = cT(:,pm0).*(rinv.^2*ones(1,length(pm0))) - 2*ep^2;
    
    cdT = zeros(Ne,M);
    cdT(:,pm) = T.r*ones(1,length(pm));
    cdT(:,pm0) = rinv*ones(1,length(pm0));
    
    T.S = cT.*T.Te(:,j-2*m+1) + cdT.*T.dTe(:,j-2*m+1);
    T.S = (T.rsc*ones(1,M)).*T.S;
    T.S(:,pm) = T.S(:,pm).*T.Pe(:,m(pm)-1 +1);
    pm00=find(m==0 & p==0);
    T.S(pr0,pm0) = 0; % First make all zero
    T.S(pr0,pm00) = ones(length(pr0),1)* ...
	  ((-1).^(j(pm00)/2+1).*(2*ep^2 + j(pm00).^2))';
  end  
      
      
      