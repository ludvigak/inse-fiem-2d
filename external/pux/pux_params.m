function PUXParams = pux_params(M, LBox, ep, P, arcL)
    
    % Percentage of parition radius when a zero-paritions will be removed
    R_ZeroPart_threshold = 0.8;
    % radius of partition ratio (1 = coverage, >1 overlap)
    R_ratio = 2.5;     
    hgrid = 2*LBox/M;

    R_part = P*hgrid;
    R = R_part/R_ratio;

    disp(['P = ' num2str(P)]);
    disp(['R = ' num2str(R)]);
    disp(['R_part = ' num2str(R_part)]);
    
    regularity = 0;
    if regularity == 0
        regularity= min(floor(sqrt(P)-0.9),5); % test
    end
    % Number of Vogel points, i.e. RBF-centres. 
    n_RBFCentres = round(min(0.8*pi/4*P^2,4*P)); 
    disp(['Number of Vogel nodes:' num2str(n_RBFCentres)])

    % Number of partitions
    n_Part = (ceil(arcL/(2*R))+1);
    
    
    
    % Pack into struct
    s = struct(...
        'R_ZeroPart_threshold', R_ZeroPart_threshold, ...
        'R_part',   R_part, ...
        'regularity', regularity,  ...
        'n_RBFCentres',         n_RBFCentres,        ...
        'n_Part',                n_Part, ...
        'M', M, ...
        'LBox', LBox, ...
        'ep', ep, ...
        'arcL', arcL );
    PUXParams = orderfields(s);
end