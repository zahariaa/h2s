function y = log2space(d1, d2, n, base)
%LOGSPACE Logarithmically spaced vector.
%   LOGSPACE(X1, X2) generates a row vector of 50 logarithmically
%   equally spaced points between octaves 2^X1 and 2^X2.  If X2
%   is pi, then the points are between 2^X1 and pi.
%
%   LOGSPACE(X1, X2, N) generates N points.
%   For N = 1, LOGSPACE returns 2^X2.
%
%   Class support for inputs X1,X2:
%      float: double, single
%
%   See also LOGSPACE, LINSPACE, COLON.
% 
% 2014-02-01 AZ Modified to be log2

%   Copyright 1984-2012 The MathWorks, Inc. 
%   $Revision: 5.11.4.6 $  $Date: 2012/02/09 20:58:06 $

switch nargin
   case 2
      n    = 50;
      base = 2;
   case 3
      base = 2;
end

if d2 == pi || d2 == single(pi) 
    d2 = log2(d2)/log2(base);
end

y = base .^ linspace(d1, d2, n);
