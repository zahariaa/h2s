function varargout = cellify(C,FREE)
% Minifunction that dumps cell contents
% 
% Code:
% C = C{:};
% 
% 2012-12-06 AZ Created

if ~exist('FREE','var') || isempty(FREE) || ~FREE
   varargout{1         }    = [C{:}];
else % FREE
   varargout(1:numel(C))    =  C    ;
end
