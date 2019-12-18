function varargout = dealvec(v)
% dealvec: Deals elements of a vector to output variables
% Ex.
% [orPref,sfPref,tfPref,orWidth,sfWidth,tfWidth] = dealvec([0 0.5 8 1 2 3]);

nv        = numel(v);
varargout = cell(1,nv);

for i = 1:nv
   varargout{i} = v(i);
end