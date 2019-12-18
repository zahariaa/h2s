function tun = tuningWrapSort(tun,wrapfactor,nbins)
% TUNINGWRAPSORT: wraps tuning bins 

if ~exist('nbins','var') || isempty(nbins)
   nbins = size(tun.bins,2);
end

if ~exist('wrapfactor','var') || isempty(wrapfactor)
   wrapfactor = [1 1 2 2];
   wrapfactor = wrapfactor(1:min(nbins,4));
elseif numel(wrapfactor)==1
   wrapfactor = repmat(wrapfactor,[nbins 1]);
end

% wrap bins
for i = 1:nbins
   tun.bins(:,i) = wrap( tun.bins(:,i), [-180 180]/wrapfactor(i) );
end

% wrap responses
for i = 1:nbins
[     tun.bins(    :,i),ix] = sort(    tun.bins(     :,  i  ));
   if nbins < size(tun.respmean,2) % 2 bins, 4 responses   
      tun.respmean(:,2*i-1) =          tun.respmean(ix,2*i-1);
      tun.respstd( :,2*i-1) =          tun.respstd( ix,2*i-1);
      tun.respmean(:,2*i  ) =          tun.respmean(ix,2*i  );
      tun.respstd( :,2*i  ) =          tun.respstd( ix,2*i  );
   else % 2/4 bins, 2/4 responses
      tun.respmean(:,  i  ) =          tun.respmean(ix,  i  );
      tun.respstd( :,  i  ) =          tun.respstd( ix,  i  );
   end
end

return
