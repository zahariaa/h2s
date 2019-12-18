function ix = sortix(varargin)
% Returns the index of sort, rather than actual sorted values
% Good for use in cellfun.
% 
% 2016-02-15 AZ Created (from minix)

switch nargin
   case 1
      [~,ix] = sort(varargin{1});
   case 2
      [~,ix] = sort(varargin{1},varargin{2}); % WILL PRODUCE ERROR
   case 3
      [~,ix] = sort(varargin{1},varargin{2},varargin{3});
end