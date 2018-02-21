function str=any2str(varargin)

% FUNCTION
%       concatenates strings and numbers into a single string. by default,
%       numbers are converted to strings with 3-digit precision. the
%       precision (number of digits) can be passed as final argument in the
%       form of a cell (to indicate that it is not to be included in the
%       concatenation.
%
% USAGE EXAMPLES
%       disp(any2str('a=',a,'b=',b,'c=',c))
%       disp(any2str('a=',a,'b=',b,'c=',c,{9}))
%
% VERSION
% updated to handle logicals by niko kriegeskorte (2016-07-30)

% unwrap if necessary
if iscell(varargin) && nargin==1
    varargin = varargin{1};
end

% if last entry is a cell, interpret as precision setting
if iscell(varargin{numel(varargin)})
    precision=cell2mat(varargin{nargin});
else
    %precision=3;
    precision='%6.3f';
end

str=[];
for argI=1:numel(varargin)
    if islogical(varargin{argI})
        bool=varargin{argI};
        bool=bool(:)';
        boolStr = repmat('0',size(bool));
        boolStr(bool)='1';
        str=[str,boolStr];
    elseif isa(varargin{argI},'char')
        str=[str,varargin{argI}];
    elseif isnumeric(varargin{argI})
        mat = varargin{argI};
        vec = mat(:);
        precision=min(numel(num2str(floor(varargin{argI}),8))+2,4);
        str=[str,num2str(vec',precision)];
    end
end