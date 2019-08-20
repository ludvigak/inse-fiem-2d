function [deg,diff,op]=RBF_QR_parse(op)
% The purpose of the function is to expand the operator information
% into the needed details
%--- op (string(:)) : Describes the operator to evaluate
%--- deg (scalar) : The degree of the operator
%--- diff (1:dim) : The number of diffs in each dimension

%--- Trim away spaces
pos = find(op~=' ');
op=op(pos);
deg=[];
diff=[];

if (op(1)=='1')
  deg=0;
elseif (op(1)=='x' | op(1)=='y' | op(1)=='z')
  deg=length(op);
  if (deg>2)
    error('Mixed derivatives of higher degree than 2 are not implemented')
  end  
  pos = find(op=='x');
  if (length(pos)>0) diff(1)=length(pos); end
  pos = find(op=='y');
  if (length(pos)>0) diff(2)=length(pos); end
  pos = find(op=='z');
  if (length(pos)>0) diff(3)=length(pos); end
elseif (op(1)=='L')
  diff=[];
  deg=2;
  if (length(op)>1)
    deg = 2*str2num(op(2:end));
  end
  if (deg>10)
    error('Hyperviscosity implemented only up to degree L^5')
  end  
end
