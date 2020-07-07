function varargout = dealvec(v)
% dealvec: Deals elements of a vector to output variables
% 2020-07-07 AZ extended to cell inputs too
% Ex.
% [orPref,sfPref,tfPref,orWidth,sfWidth,tfWidth] = dealvec([0 0.5 8 1 2 3]);
% [a,b,c] = dealvec({[3 5] 'a' 100});

nv = numel(v);

if iscell(v)
   varargout = v(1:min(nargout,nv));
else
   varargout = cell(1,nargout);
   for i = 1:min(nargout,nv)
      varargout{i} = v(i);
   end
end
