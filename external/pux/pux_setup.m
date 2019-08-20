function PUXData = pux_setup(xe, inside, z_grid, z_equi, params)
% Assumptions that caller must take care of:
% - no zero partitions in concave areas    
    
% Set PUX-parameters
R_ZeroPart_threshold = params.R_ZeroPart_threshold;
R_part = params.R_part;
regularity = params.regularity;
n_RBFCentres = params.n_RBFCentres;
n_Part = params.n_Part;
M = params.M;
LBox = params.LBox;
ep = params.ep;

%% 
% Create PUX grid on box D



% Evalutaion points
% Number of evaluation points (not sure I'll need it.)
n_Boxpnts = size(xe,1);
% Create a NeighborSearcher object for K-nearest neighbors search.
NeighborSearcher = createns(xe);
%% List nodes that are inside or outside Omega.
% The boundary of Omega is given by the equally spaced points (w.r.t.
% arclength) x_bdry.

% Mark the points in xe as inside or outside Omega. If point xe(i,:) is
% inside element then idx_xe_in_boolean(i) = 1,m oterhwise 0.
%idx_xeInOmega_boolean = inpolygon(xe(:,1), xe(:,2), real(z), imag(z));
idx_xeInOmega_boolean = inside;
% List containing the indices of points xe that are outside Omega.
idx_xeOutOmega = find(~idx_xeInOmega_boolean);

%inner_crs = i_crs(idx_in_logic(i_crs));

%% Setup partitions
[zD,~] = fft_diff(z_equi ,[1,numel(z_equi)]);
n_equi = -1i*zD./abs(zD);
[idx_InterpPartCentres, idx_ZeroPartCentres, R_ZeroPart] = ...
    pux_partitions(xe, inside, z_equi, n_equi, z_grid, M, LBox, R_part, R_ZeroPart_threshold);

ZeroPartCentres = xe(idx_ZeroPartCentres,:);
n_ZeroPart = size(ZeroPartCentres,1);
InterpPart = xe(idx_InterpPartCentres,:);
n_InterpPart = n_Part;
R_InterpPart = R_part;
% Total number of partitions
n_Part = n_ZeroPart + n_InterpPart;
% Locate index in xe for points withing R_InterpPart of InterpPart.
% Sufficient to do this for one interpolation partition and then shift.
[cell_idx_xeInInterpPart_1,Dist_xe2InterpPartCentre_1] = rangesearch(NeighborSearcher,InterpPart(1,:),R_InterpPart);
[idx_xeInInterpPart_stencil, idx_sort_ascending_order] = sort(cell_idx_xeInInterpPart_1{1}-idx_InterpPartCentres(1));
Dist_xe2InterpPartCentre_stencil = Dist_xe2InterpPartCentre_1{1};
Dist_xe2InterpPartCentre_stencil = Dist_xe2InterpPartCentre_stencil(idx_sort_ascending_order)';
n_pntsInInterpPart = length(idx_xeInInterpPart_stencil);
% Create Vogel-points. These will be used as RBF-centers.
vogel_param = pi*(3-sqrt(5));
RBFCentres = R_InterpPart * repmat(sqrt(1:n_RBFCentres)'/sqrt(n_RBFCentres),1,2).* [cos((1:n_RBFCentres)'*vogel_param),sin((1:n_RBFCentres)'*vogel_param)];

n_extensionPnts = round(sum(n_InterpPart) * n_pntsInInterpPart*2/3);


%% RBF--Direct 
% xe1 = xe(idx_inIp + idx_Ip(1),:);
%     xv1 = x_rbf+repmat(Ip(1,:),n_v,1);
%     
%     [x1,x2] = ndgrid(xe1(:,1),xv1(:,1));
%     [y1,y2] = ndgrid(xe1(:,2),xv1(:,2));
%     r2 = (x1-x2).^2 + (y1-y2).^2;
%     A = exp(-params.ep^2*r2);
%% RBF--QR
% Create PUX-matrix for the stencil interpolation partition.
A_PUX = RBF_QR_diffmat_2D('1', xe(idx_xeInInterpPart_stencil + idx_InterpPartCentres(1),:), [RBFCentres(:,1)+InterpPart(1,1), RBFCentres(:,2)+InterpPart(1,2)], ep, 2);



% Bad measure
disp(['Overdetermined with factor ' num2str(size(A_PUX,1)/size(A_PUX,2))])

%% Choose RBF to use as building block for weight function.

switch regularity
    case 1
        disp('Shepard with Wu1')
        Phi = @(r) 1/2*(1-r).^2.*(2+r); % Wu C0(R3)
    case 2
        disp('Shepard with Wu2')
        Phi = @(r) 1/8*(1-r).^3.*(8+9*r+3*r.^2); % Wu C0(R5)
    case 3
        disp('Shepard with Wu3')
%       Phi = @(r) 1/16*(1-r).^4.*(16+29*r+20*r.^2+5*r.^3); % Wu C0(R7)
        Phi = @(r) 1/4*(1-r).^4.*(4+16*r+12*r.^2+3*r.^3); % Wu C2(R3)
    case 4
        disp('Shepard with Wu4')
        Phi = @(r) 1/8*(1-r).^5.*(8+40*r+48*r.^2+25*r.^3+5*r.^4); % Wu C2(R5)
    case 5
        disp('Shepard with Wu5')
        Phi = @(r) 1/6*(1-r).^6.*(6+36*r+82*r.^2+72*r.^3+30*r.^4+5*r.^5); % Wu C4(R3)
    otherwise
        disp('Shepard with Wendland C2')
        Phi = @(r) (1-r).^4.*(4*r+1); % Wendland C2(R3)
end

%% Evaluate an RBF associated with an interpolation partition once, then reuse.
RBF_evalInInterpPart_stencil = Phi(Dist_xe2InterpPartCentre_stencil/R_InterpPart);

%% Loop over interpolation partitions

n_RBFevals = n_Part*n_pntsInInterpPart;


% Vectors for data used for building extension. The names I, J and S is
% based on the input for MATLABS sparse:
% S = sparse(i,j,s,m,n) uses vectors i, j, and s to generate an
%     m-by-n sparse matrix such that S(i(k),j(k)) = s(k)

% If, Jf and Sf stores information about the local extensions.
If = zeros(n_extensionPnts,1); 
Jf = zeros(n_extensionPnts,1);
% Iw, Jw and Sw stores information about the local weights.
Iw = zeros(n_RBFevals,1); Jw = zeros(n_RBFevals,1); Sw = zeros(n_RBFevals,1);

% Counters used to fill in If, Jf, Sf, Iw, Jw and Sw.
idx_endf = 0;
idx_endw = 0;
msgid = 'MATLAB:nearlySingularMatrix';
warning('off',msgid);
msgid2 = 'MATLAB:rankDeficientMatrix';
warning('off',msgid2);


minj = 10^16; % Stores the number of interpolation points for the
% interpolation partitions with the fewest. Used to determine if the least
% squares systems are sufficiently "overdetermined".
maxj = 0; % Same as minj, but stores the maximum amount of interpolation
% points. Currently not used.

for i = 1:n_InterpPart
    % Global indices
    idx_xeInInterpPart_i = idx_InterpPartCentres(i) + idx_xeInInterpPart_stencil; % Indices of point in xe within R_Ip of parition center Ip(i).
    idx_xeInOmegaInInterpPart_i = idx_xeInInterpPart_i(idx_xeInOmega_boolean(idx_xeInInterpPart_i));
    idx_xeInInterpPartOutOmega = idx_xeInInterpPart_i(~idx_xeInOmega_boolean(idx_xeInInterpPart_i));
    n_local_evals = length(idx_xeInInterpPartOutOmega);
    
    maxj = max(length(idx_xeInOmegaInInterpPart_i),maxj);
    minj = min(length(idx_xeInOmegaInInterpPart_i),minj);

    If(idx_endf+1:idx_endf+n_local_evals) = idx_xeInInterpPartOutOmega;
    Jf(idx_endf+1:idx_endf+n_local_evals) = i;

    % Save contribution to weight function
    Iw(idx_endw+1:idx_endw+n_pntsInInterpPart) = idx_xeInInterpPart_i;
    Jw(idx_endw+1:idx_endw+n_pntsInInterpPart) = i;
    Sw(idx_endw+1:idx_endw+n_pntsInInterpPart) = RBF_evalInInterpPart_stencil;
    idx_endf = idx_endf+n_local_evals;
    idx_endw = idx_endw+n_pntsInInterpPart;
    
    
end

% cell_idx_xeOutOmegaInInterpPart is obtained by range search from 
%partition centers on uniform grid. Thus they are not the same amount as 
%idx_local_xeOutOmega_boolean. idx_local_xeOutOmega_boolean
% comes from inpolygon. 
If = If(1:idx_endf); 
Jf = Jf(1:idx_endf); 

idx_xeOutOmegaInInterpPart = sort(unique(If));


warning('on',msgid)
warning('on',msgid2) 
 beta_min = minj/n_RBFCentres;
if beta_min<2
    warning(['Least-squares oversampling less than two (' num2str(beta_min) ')!']);
end
disp(['beta_min: ' num2str(beta_min)])
%% Calculate contribution to weight function for zero-partitions
for i = 1:n_ZeroPart    
    idx_xeInZeroPart_i = idx_ZeroPartCentres(i) + idx_xeInInterpPart_stencil;
    if  max(idx_xeInZeroPart_i)> n_Boxpnts || min(idx_xeInZeroPart_i)< 0
        idx_insideB = intersect(find(idx_xeInZeroPart_i < n_Boxpnts+1),find(idx_xeInZeroPart_i > 0));
        n_xeInZeropart_i = length(idx_insideB);
        Dist_xe2ZeroPartCentre_i = Dist_xe2InterpPartCentre_stencil(idx_insideB);
        Iw(idx_endw+1:idx_endw+n_xeInZeropart_i) = idx_xeInZeroPart_i(idx_insideB);
        Jw(idx_endw+1:idx_endw+n_xeInZeropart_i) = i+n_InterpPart;
        Sw(idx_endw+1:idx_endw+n_xeInZeropart_i) = Phi(Dist_xe2ZeroPartCentre_i/R_ZeroPart(i));
        idx_endw = idx_endw+n_xeInZeropart_i;
    elseif norm(R_ZeroPart(i)- R_InterpPart) > 1e-14
      
        Dist_xe2ZeroPartCentre_i = find(Dist_xe2InterpPartCentre_stencil<=R_ZeroPart(i));
        idx_xeInZeroPart_i = idx_xeInZeroPart_i(Dist_xe2ZeroPartCentre_i);
        Dist_xe2ZeroPartCentre_i = Dist_xe2InterpPartCentre_stencil(Dist_xe2ZeroPartCentre_i);
        n_xeInZeropart_i = length(Dist_xe2ZeroPartCentre_i);
        Iw(idx_endw+1:idx_endw+n_xeInZeropart_i) = idx_xeInZeroPart_i;
        Jw(idx_endw+1:idx_endw+n_xeInZeropart_i) = i+n_InterpPart;
        Sw(idx_endw+1:idx_endw+n_xeInZeropart_i) = Phi(Dist_xe2ZeroPartCentre_i/R_ZeroPart(i));
        idx_endw = idx_endw+n_xeInZeropart_i;
    else
        Iw(idx_endw+1:idx_endw+n_pntsInInterpPart) = idx_ZeroPartCentres(i) + idx_xeInInterpPart_stencil;
        Jw(idx_endw+1:idx_endw+n_pntsInInterpPart) = i+n_InterpPart;
        Sw(idx_endw+1:idx_endw+n_pntsInInterpPart) = RBF_evalInInterpPart_stencil;
        idx_endw = idx_endw+n_pntsInInterpPart;
    end
    
end


%% Build weight matrix W

% Remove zero-valued elements, these occur because all paritions do not
% hold n_pntsInInterpPart points.
Iw = Iw(1:idx_endw); Jw = Jw(1:idx_endw); Sw = Sw(1:idx_endw);
W = sparse(Iw,Jw,Sw,n_Boxpnts,n_Part);
W = W(idx_xeOutOmegaInInterpPart,:);
sumW = sum(W,2);
W = spdiags(1./sumW,0,length(idx_xeOutOmegaInInterpPart),length(idx_xeOutOmegaInInterpPart))*W;

function s=makestruct(idx_xeInOmega_boolean, Ip, Op, R_Ip, R_Op,A_PUX,W_PUX,If,Jf,n_Boxpnts,n_Part,n_InterpPart,n_extensionPnts,idx_xeOutOmegaInInterpPart,idx_InterpPartCentres,idx_xeInInterpPart_stencil)
    s = struct('A',A_PUX,'W',W_PUX,'If',If,'Jf',Jf,'Ip',Ip,'Op', ...
               Op,'R_Ip',R_Ip,'R_Op',R_Op,'n_Boxpnts',n_Boxpnts,'n_Part', ...
               n_Part,'n_InterpPart',n_InterpPart,'n_extensionPnts', ...
               n_extensionPnts, 'idx_xeOutOmegaInInterpPart', ...
               idx_xeOutOmegaInInterpPart,'idx_InterpPartCentres', ...
               idx_InterpPartCentres, 'idx_xeInInterpPart_stencil', ...
               idx_xeInInterpPart_stencil,'idx_xeInOmega_boolean', ...
               idx_xeInOmega_boolean);
end

PUXData = makestruct(idx_xeInOmega_boolean, InterpPart, ZeroPartCentres, R_InterpPart, R_ZeroPart,A_PUX,W,If,Jf,n_Boxpnts,n_Part,n_InterpPart,n_extensionPnts,idx_xeOutOmegaInInterpPart,idx_InterpPartCentres,idx_xeInInterpPart_stencil);

end
