% tolD, Check how to use that now
% Antal vektorer i Q och de Qnew som l??ggs till har r??tt storlek. Q.q ??r
% alltid n??st plats d??r det ska l??ggas in.

% Fixa order, cols, qused in more ways than one
function [Q,R] = IncQR(newC,bf,bj,tol,Q,R)
% This function performs an incremental QR-factorization
% with a block pivoting strategy specially designed for the RBF-QR method.
%--- newC(1:N,:) : The new columns to incorporate into the QR-fac
%--- bf         : The scaling factor between two blocks  
%--- bj(1:)     : The block number of each column. Columns should come in
%                 whole blocks.
%--- tol        : Determines when to cut out columns   
%--- Q (struct) : Information about Q and pivoting
%--- R (struct) : The resulting columns in R

mp=eps; %Machine precision
zz = 10*mp;
N = size(newC,1);
%
% Check if first time operating on C
%
if (nargin<=4) % First time
  R.cn = 0; % The number of handled/seen columns
  R.def_j = [];
  R.def_col =zeros(N,0);
  R.def_R = zeros(N,0);
  R.R = zeros(N,0);
  R.order = [];
  R.def_ord = [];
  R.def_pow = [];
  Q.q = 1; % The next q-vector
  Q.Q = zeros(N,0);
end

cn = size(newC,2);
cols = R.cn + (1:cn); % Original column numbers
R.cn = R.cn + cn;     % New starting point for next time

% --- Compute upper R-values for all new columns
R0 = Q.Q'*newC;

% --- If we have not yet filled Q 
if (Q.q <= N) 
  
  % --- Remove the Q-component from all columns
  newC = newC - Q.Q*(Q.Q'*newC);
  
  % --- Work block by block through the matrix (assuming ordered)
  for j=bj(1):bj(end)
    R.def_num(max(1,j)) = 0; % Never defer for j=0
    % --- Locate the current block
    pos = find(bj==j);
    n = length(pos);
    % --- Check for deferred columns
    dpos = find(R.def_j==j);
    nb = n + length(dpos);
    % --- Store actual column numbers for later use
    order(1:n) = cols(pos);
    order(n+1:nb) = R.def_ord(dpos);
    % --- QR-factorize the block with pivoting
    [Qnew,Rnew,Enew]=qr([newC(:,pos) R.def_col(:,dpos)],0);
    % --- Check the newly computed diagonal elements for significance
    dR = abs(diag(Rnew));
    % --- Plot the latest R_ii elements against block index
    %figure(10),H=semilogy(j,dR,'r+'); hold on
    %set(H,'LineWidth',2)
    R.origR{j+1}=dR;
    % --- Look for elements that drop significantly in magnitude
    mag = log10(dR);
    ll=length(mag);
    if (ll<nb)
      mag(ll+1:nb) = mag(ll); % We have passed column N
    end
    diff = mag(1:end-1)-mag(2:end);
    pp = find(diff > tol); % As default tolerance, I suggest 2 (100 times)
    pp = [pp;nb]; % Later used as intervals, adding end 
    % --- Collect the part of R corresponding to old q-vectors
    upper = [R0(:,pos) R.def_R(1:Q.q-1,dpos)];

    % --- Split into two parts. The one to keep and the one to defer  
    upper = upper(:,Enew(1:pp(1))); % Sort and cut
  
    % ---  Remove the parts that will not be used
    Qnew = Qnew(:,1:min(N,pp(1)));
    Rnew = Rnew(1:min(N,pp(1)),1:pp(1));
    if (Q.q-1+pp(1) > N) % We have more columns than we need left
      rows = N-Q.q+1;
      Rnew = Rnew(1:rows,:);
      Qnew = Qnew(:,1:rows);
    end  
    ddR = diag(Rnew);
    if (size(ddR,1)==size(ddR,2))
      ddR=ddR(1); % Special ifjust one row
    end
    %figure(10),H=semilogy(j,abs(ddR),'bo'); hold on
    %set(H,'MarkerSize',14,'LineWidth',2)
    R.selectedR{j+1} = abs(ddR); 
    % --- Reorthogonalise the new basis vectors (otherwise drift kills accuracy)
    Qnew = Qnew - Q.Q*(Q.Q'*Qnew);
    %--- Renormalize also?
    qnorm=sqrt(diag(Qnew'*Qnew));
    Qnew = Qnew*diag(1./qnorm);
    % --- Update Q and R
    R.R = [R.R [upper; Rnew; zeros(N-Q.q+1-size(Rnew,1),pp(1))]];
    Q.Q = [Q.Q Qnew];
    Q.q = Q.q + size(Qnew,2);

    R.order = [R.order order(Enew(1:pp(1)))];
    
    % --- Apply the new q-vectors to the remaining vectors
    R0 = [R0; [zeros(size(Qnew,2),pos(end)) Qnew'*newC(:,pos(end)+1:end)]];
    newC(:,pos(end)+1:end) = newC(:,pos(end)+1:end) - Qnew*(Qnew'*newC(:,pos(end)+1:end));
    % --- Handle the deferred columns. Assuming tol corresponds to blockdiff 
    if (length(pp)>1)       % There are columns to defer
      % --- diff = -log10(bf^p) = - p log10(bf)
      if (bf<1)
        p = -diff(pp(1:end-1))/log10(bf); % Default bf=ep^2
      else
        p = ones(size(diff(pp(1:end-1)))); % We can move anything if ep is large
      end  
      pow = zeros(1,nb-pp(1));
      for ll=1:length(pp)-1
        pow(pp(ll)+1:pp(ll+1)) = floor(sum(p(1:ll)));
      end
      % --- First handle all new deferred columns
      p1 = find(Enew(pp(1)+1:end) <= n); p1 = p1 + pp(1);
      loc = Enew(p1);      ll = length(loc);
      R.def_num(j)=R.def_num(j)+ll; % Counter
      R.def_col = [R.def_col newC(:,pos(loc))];
      R.def_R = [R.def_R [R0(:,pos(loc)); zeros(N-Q.q+1,ll)]];
      R.def_ord = [R.def_ord cols(pos(loc))];
      R.def_pow = [R.def_pow pow(p1)];
      if (mag(p1)<=log10(10*zz))
%         'Zero'
        pow(p1) = pow(p1) + 5; % Zero, move forward
      end
      if (Q.q > N)
%         'Already at the end'
        pow(p1)=1; % Take these into the final R at the end
      end
      R.def_j = [R.def_j j+pow(p1)];
      % --- Go through old columns to redefer
      p2 = find(Enew(pp(1)+1:end) > n); p2 = p2 + pp(1);
      for k=1:length(p2) 
        loc = Enew(p2(k))-n;
        % --- Check if we can move it even further. It is small also
        % here.
%         'Redefer?'
        pow(p2(k)) = pow(p2(k)) - R.def_pow(dpos(loc));
        if (mag(p2(k)) <= log10(10*zz))
%           'Zero'
          pow(p2(k)) = 5; % Zero, move forward
        end
        if (pow(p2(k))>0)
%           'Moving'
          R.def_pow(dpos(loc)) = R.def_pow(dpos(loc)) + pow(p2(k));
          R.def_j(dpos(loc)) = j+pow(p2(k));
        else
          R.def_pow(dpos(loc));
          powerval = pow(p2(k));
          magval = mag(p2(k));
          %pause
          % MAY WANT TO CHANGE THIS BACK TO AN ERROR
          warning(['Found an undeferable column. This case is not defined ' ...
                 'yet'])
          % Probably we should then decide to add the column
          pow(p2(k))=1;
          R.def_pow(dpos(loc)) = R.def_pow(dpos(loc)) + pow(p2(k));
          R.def_j(dpos(loc)) = j+pow(p2(k));
        end        
      end
    end
    % --- Apply the new Q to the deferred columns
    nq = size(Qnew,2);
    R.def_R(Q.q-nq:Q.q-1,:) = Qnew'*R.def_col;
    R.def_col = R.def_col - Qnew*(Qnew'*R.def_col);
    
    
    if (size(R.R,2) >=N)
      % --- Put the deferred columns back into R at the end.
      pos = find(R.def_j > j);
      % --- Make nearly zero into exactly zero
      for q=1:length(pos)
        where = find(abs(R.def_R(:,pos(q)))<=10*zz);
        R.def_R(where,pos(q)) = 0;
      end	
      R.R = [R.R R.def_R(:,pos)]; % possibly sort them
      R.order = [R.order R.def_ord(pos)];
      % --- Put in the columns that were not treated yet.
      pos = find(bj >j);
      rz = size(R0,1);
      R.R = [R.R [R0(:,pos); Q.Q(:,rz+1:end)'*newC(:,pos)]];
      R.order = [R.order cols(pos)];
      break
    end
  end
else
  % --- If we have filled all columns of Q, just project and join
  R.R = [R.R R0];
  R.order = [R.order cols];
end
origR=R.origR;
selectedR= R.selectedR;

