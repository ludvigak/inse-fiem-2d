function [A,T]=RBF_QR_mat_2D(Psi,op,var)
% The purpose of this function is to return an evaluation matrix A
% for the given operator op. The elements A_{ij} are given by
% op(psi_j(xe(i))), where psi is the RBF-QR basis and xe are the
% evaluation points.
%--- Psi (struct) : Defines the basis functions
%--- op (string) : Defines the operator
%--- var (double(1:N,1:2) or struct) : Either the evaluation points
%---    xe or a structure containing precomputed data including xe
%
%--- Call patterns
%--- [A,T]=RBF_QR_mat_2D(Psi,op,xe) First time for these points (xe)
%--- [A,T]=RBF_QR_mat_2D(Psi,op,T)  Subsequent times with same xe
%--- [A]=RBF_QR_mat_2D(Psi,op,xe) Just one call with nice interface

%--- Find out which call sign we got

if isstruct(var)
    T=var;
    xe = T.xe;
else
    T=[];
    xe=var;
end

[deg,diff,op] = RBF_QR_parse(op);
%--- Precompute functions that are needed for evaluation of T_{j,m}
ep = Psi.ep;
Rt = Psi.Rt;
N = size(Psi.Rt,1);
order = Psi.columns;
M = length(order);
j = Psi.j(order);
m = Psi.m(order);
cs = Psi.cs(order);
p = Psi.p(order);
T = RBF_QR_precomp_2D(j,m,p,ep,xe,deg,T);

%--- Evaluate the basis functions
if (deg==0)
    V = (T.rsc*ones(1,M)).*T.Pe(:,m+1).*T.Te(:,j-2*m+1);
    pos = find(cs == 1);    V(:,pos) = V(:,pos).*T.Hec(:,2*m(pos)+p(pos));
    pos = find(cs == -1);   V(:,pos) = V(:,pos).*T.Hes(:,2*m(pos)+p(pos));
    
    A = V(:,1:N) + V(:,N+1:end)*Rt.';
elseif (deg==1)
    co1 = [T.Hec(:,1) T.Hes(:,1)];
    co2 = [1 -1];
    if (op=='x')
        i1 = 1; i2 = 2;
    elseif (op=='y')
        i1 = 2; i2 = 1;
    end
    A = (co1(:,i1)*ones(1,M)).*T.G;
    V = (co1(:,i2)*ones(1,M)).*T.F;
    pos = find(cs == 1);
    A(:,pos) = A(:,pos).*T.Hec(:,2*m(pos)+p(pos)) + ...
        co2(i1)*V(:,pos).*T.Hes(:,2*m(pos)+p(pos));
    pos = find(cs == -1);
    A(:,pos) = A(:,pos).*T.Hes(:,2*m(pos)+p(pos)) + ...
        co2(i2)*V(:,pos).*T.Hec(:,2*m(pos)+p(pos));
    
    A = A(:,1:N) + A(:,N+1:end)*Rt.';
elseif (deg==2)
    if (op(1)=='L')
        A = T.Q+T.S;
        pos = find(cs == 1);
        A(:,pos) = A(:,pos).*T.Hec(:,2*m(pos)+p(pos));
        pos = find(cs == -1);
        A(:,pos) = A(:,pos).*T.Hes(:,2*m(pos)+p(pos));
    else
        % g(theta)=cos((2m+p)theta) => h(theta)=-sin()
        % g(theta)=sin((2m+p)theta) => h(theta)= cos()
        co1 = [T.Hec(:,1).^2 T.Hes(:,1).^2 T.Hec(:,1).*T.Hes(:,1) T.Hec(:,2)];
        co2 = [1 1 -1];
        co3 = [2 -2 -1];
        if (op=='xx')
            i1 = 1;    i2 = 2;    i3 = 3;
        elseif (op=='xy')
            i1 = 3;    i2 = 3;    i3 = 4;
        elseif (op=='yy')
            i1 = 2;    i2 = 1;    i3 = 3;
        end
        A = (co1(:,i1)*ones(1,M)).*T.Q + (co2(i1)*co1(:,i2)*ones(1,M)).*T.S;
        V = (co3(i1)*co1(:,i3)*ones(1,M)).*T.R;
        
        pos = find(cs == 1);
        A(:,pos) = A(:,pos).*T.Hec(:,2*m(pos)+p(pos)) - ...
            V(:,pos).*T.Hes(:,2*m(pos)+p(pos));
        pos = find(cs == -1);
        A(:,pos) = A(:,pos).*T.Hes(:,2*m(pos)+p(pos)) + ...
            V(:,pos).*T.Hec(:,2*m(pos)+p(pos));
    end
    
    A = A(:,1:N) + A(:,N+1:end)*Rt.';
end
