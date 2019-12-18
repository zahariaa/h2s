function data = wrap(data,xrange)
% wrap.m: wraps circular data (COLUMNWISE) to keep it in a specified range
% e.g.,
% wrap(data,[0 360])
% wrap(data,[-pi pi])
% 
% 2014-09-03 AZ Created

if diff(xrange)==0 || abs(diff(xrange))==Inf;   return;   end

% Deal with row vectors, matrices
FLIP = false;
if all(size(data) > 1)
   % recurse on columns
   for i = 1:size(data,2)
      data(:,i) = wrap(data(:,i),xrange);
   end
   return
elseif size(data,1)==1 && size(data,2) > 0 && mod(size(data,2),1)==0 %==isrow(data)
   FLIP = true;
   data = data';
end

% amount of shift needed
bump = diff(xrange);

% keep positions of extrema
preserve = [data==xrange(1) data==xrange(2)];

% apply shift
data(data < xrange(1)) = data(data < xrange(1)) + bump;
data(data > xrange(2)) = data(data > xrange(2)) - bump;

% ugly heuristic for forcing sequential bounds
if any(sum(preserve)==2)
   if     sign(sum(diff(sign(diff(data))))) == 1
      data(  1) = xrange(1);
   elseif sign(sum(diff(sign(diff(data))))) == -1
      data(end) = xrange(2);
   end
end

% restore bounds if duplicate
if     any(preserve(:,1))
   data(sum([data==xrange(1) preserve(:,1)],2)==1) = xrange(2);
elseif any(preserve(:,2))
   data(sum([data==xrange(2) preserve(:,2)],2)==1) = xrange(1);
end

% Recurse if job is not done!
if any(data < xrange(1)) || any(data > xrange(2))
   data = wrap(data,xrange);
end

if FLIP;  data = data';  end

end