function slop = reshape_sloppy(x,newsz)
% slop = reshape_sloppy(x,newsz)
% Reshapes x to an imperfect size padded with NaNs
%
% 2018-04-24 AZ Created

%% Preliminaries
oldsz = size(x);
if numel(newsz)==1
   newsz(2) = ceil(prod(oldsz)/newsz);
end

%% RESHAPE!
newrows = prod(oldsz)-prod(newsz);
if newrows>0
   newsz(1) = newsz(1) + ceil(newrows/oldsz(1));
end
slop = NaN(newsz);
slop(1:prod(oldsz)) = x';
slop = slop';

return

