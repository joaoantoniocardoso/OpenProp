% -------------------------------------------------------------------------
%  pchipslopes - First derivative slopes for shape-preserving Hermite cubic
%
%  pchipslopes(x,v) computes d(i) = dvdx(x(i)).
%
%  Slopes at interior points
%       delta = diff(v)./diff(x)
%
%       d(i)  = 0 if delta(i-1) and delta(i) have opposites signs or either is zero.
%       d(i)  = weighted harmonic mean of delta(i-1) and delta(i) if they have the same sign.
%
% ASSUMES v = v(x) is 1D data

function d = pchipslopes(x,v)
   h     = diff(x);
   delta = diff(v)./h;
   
   Nx   = length(h)+1;
   d    = zeros(size(h));
   i    = find(sign(delta(1:Nx-2)).*sign(delta(2:Nx-1))>0)+1;
   w1   = 2*h(i)+h(i-1);
   w2   = h(i)+2*h(i-1);
   d(i) = (w1+w2)./(w1./delta(i-1) + w2./delta(i));

%  Slopes at endpoints

   d(1)  = pchipend(h(1),h(2),delta(1),delta(2));
   d(Nx) = pchipend(h(Nx-1),h(Nx-2),delta(Nx-1),delta(Nx-2));

% -------------------------------------------------------------------------

function d = pchipend(h1,h2,del1,del2)
%  Noncentered, shape-preserving, three-point formula.
   d = ((2*h1+h2)*del1 - h1*del2)/(h1+h2);
   if sign(d) ~= sign(del1)
      d = 0;
   elseif (sign(del1)~=sign(del2))&(abs(d)>abs(3*del1))
      d = 3*del1;
   end
