function varargout = maxix(varargin)
% Returns the index of max, rather than actual max value
% Good for use in cellfun.
% 
% 2013-05-31 AZ Created (from minix)

switch nargout
   case 2
      [varargout{1},varargout{2}] = find(varargin{1}==max(varargin{1}(:)));
      return
   otherwise % do nothing
end

switch nargin
   case 1
      [~,ix] = max(varargin{1});
   case 2
      [~,ix] = max(varargin{1},varargin{2}); % WILL PRODUCE ERROR
   case 3
      [~,ix] = max(varargin{1},varargin{2},varargin{3});
end

varargout{1} = ix;
return

