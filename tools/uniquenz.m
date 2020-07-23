function varargout = uniquenz(varargin)
% uniquenz: Minifunction that is the same as unique, but operates only on 
%    nonzero entries.
% 
% 2020-07-21 AZ Created from uniquenan

switch nargout
   case 3
      [varargout{1},varargout{2},varargout{3}] = unique(varargin{1}(~~varargin{1}),varargin{2:end});
   case 2
      [varargout{1},varargout{2}] = unique(varargin{1}(~~varargin{1}),varargin{2:end});
   otherwise
      varargout{1} = unique(varargin{1}(~~varargin{1}),varargin{2:end});
end

return
