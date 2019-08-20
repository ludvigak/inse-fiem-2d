function [idx_InterpPartCentres, idx_ZeroPartCentres, R_ZeroPart] = ...
    pux_partitions(xe, inside, z_equi, n_equi, z_grid, M, LBox, R_part, R_ZeroPart_threshold)

% Original code for meshgrid, I think I have changed in the important places
gridtype = 'ndgrid'; 

idx_xeInOmega_boolean = inside;
z = z_grid;
ZeroPartCentres_nonuniform = z_equi;
ZeroPartCentres_nonuniform_normal = n_equi;
n_Part = numel(z_equi)/2;

% [ZeroPartCentres_nonuniform,~,~] = Traparcldisc2(curve,n_Part*2); % Create new uniformly distributed points w.r.t. arclength.
% %D_ZeroPartCentres_nonuniform contains the values of the derivative of the
% %paramterization of the curve evaluated at the 2*n_Part points.
% [D_ZeroPartCentres_nonuniform,~] = fft_diff(ZeroPartCentres_nonuniform ,[1,n_ZeroPart]);
% ZeroPartCentres_nonuniform_normal = -1i*D_ZeroPartCentres_nonuniform./abs(D_ZeroPartCentres_nonuniform);

%% Find centers for interpolation partitions InterPart
% We find the centres for the n_p interpolation partitions InterpPart by
% sampling uniformly with respect to arclength along the boundary. For each
% such point the closest uniform grid point inside Omega is found and set
% to be a partition centre. Since 2*n_p zero-partitions centres are used we
% can sample them at the same time.

%h = abs(xe(2,2)-xe(1,2));
h = norm(xe(1,:)-xe(2,:));

% clf()
% plot(real(z_grid), imag(z_grid), '-k')
% hold on
% plot(xe(inside, 1), xe(inside, 2), '.b')
% plot(xe(~inside, 1), xe(~inside, 2), '.r')

InterpPartCentres_nonuniform = [real(ZeroPartCentres_nonuniform(1:2:end)),imag(ZeroPartCentres_nonuniform(1:2:end))];
idx_InterpPartCentres = zeros(n_Part,1);
for i = 1:n_Part
    % off-grid center location 
    cx = InterpPartCentres_nonuniform(i,1);
    cy = InterpPartCentres_nonuniform(i,2);
    % fuzzy indices
    
    % Using meshgrid
    switch gridtype
      case 'meshgrid'
        ifx = (cx+LBox)/h;
        ify = (cy+LBox)/h;        
        idx_InterpPartCentres_uniform = ...
            [floor(ifx)*M + floor(ify)+1;
             floor(ifx)*M + ceil(ify)+1;
             ceil(ifx)*M + floor(ify)+1;
             ceil(ifx)*M + ceil(ify)+1
            ];
      case 'ndgrid'
        ifx = 1 + (cx+LBox)/h;
        ify = 1 + (cy+LBox)/h;
        siz = [M, M];
        idx_InterpPartCentres_uniform = ...
            [sub2ind(siz, floor(ifx), floor(ify))
             sub2ind(siz, floor(ifx), ceil(ify))
             sub2ind(siz, ceil(ifx), floor(ify))
             sub2ind(siz, ceil(ifx), ceil(ify))
                   ];
        
    end       
    
    % plot(cx, cy, 'xb')
    % cand = xe(idx_InterpPartCentres_uniform, :)
    % plot(cand(:, 1), cand(:, 2), 'or')
    
    idx_InterpPartCentres_candidates = ...
        find(idx_xeInOmega_boolean(idx_InterpPartCentres_uniform)); 
    assert(~isempty(idx_InterpPartCentres_candidates), 'No interior candidates!')
    %Requires that the box E is symmetric around origo.
    idx_InterpPartCentres(i) = ...
        idx_InterpPartCentres_uniform(idx_InterpPartCentres_candidates(1));
end

R_InterpPart = R_part;
disp(['Rp = ' num2str(R_InterpPart)]);
%% Create zero-partitions Op

n_ZeroPart = n_Part*2;

ZeroPartCentres_nonuniform = [real(ZeroPartCentres_nonuniform),imag(ZeroPartCentres_nonuniform)];
ZeroPartCentres_nonuniform = [ZeroPartCentres_nonuniform(:,1) + (R_InterpPart + h/2)*real(ZeroPartCentres_nonuniform_normal), ZeroPartCentres_nonuniform(:,2) + (R_InterpPart + h/2)*imag(ZeroPartCentres_nonuniform_normal)];


idx_ZeroPartCentres = round((ZeroPartCentres_nonuniform+LBox)/h);
switch gridtype
  case 'meshgrid'
    idx_ZeroPartCentres = idx_ZeroPartCentres(:,1) * M + idx_ZeroPartCentres(:,2) + 1;
  case 'ndgrid'
    idx_ZeroPartCentres = idx_ZeroPartCentres(:,2) * M + idx_ZeroPartCentres(:,1) + 1;
end
ZeroPartCentres = xe(idx_ZeroPartCentres,:);

% Some partition centres might end inside Omega, happens in concave ares. These are removed.
%idx_ZeroPartCentresInOmega_boolean = inpolygon(ZeroPartCentres(:,1),ZeroPartCentres(:,2),real(z),imag(z));
idx_ZeroPartCentresInOmega_boolean = idx_xeInOmega_boolean( idx_ZeroPartCentres );
ZeroPartCentres = ZeroPartCentres(~idx_ZeroPartCentresInOmega_boolean,:);
idx_ZeroPartCentres = idx_ZeroPartCentres(~idx_ZeroPartCentresInOmega_boolean);

% Set radius for zero partitions, remove those that are too small.
[~,Dist_bdry2ZeroPartCentres] = knnsearch([real(z) imag(z)],ZeroPartCentres);
idx_ZeroPart_keep = find(Dist_bdry2ZeroPartCentres>R_InterpPart*R_ZeroPart_threshold);
R_ZeroPart = Dist_bdry2ZeroPartCentres(idx_ZeroPart_keep);
idx_ZeroPartCentres = idx_ZeroPartCentres(idx_ZeroPart_keep);
% If a zero-partition has greater than R_Ip radius, set its radius to R_Ip.

R_ZeroPart(R_ZeroPart>R_InterpPart) = R_InterpPart;
