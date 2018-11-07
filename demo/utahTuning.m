function tun = utahTuning(data)

bins     = vectify(unique(data.stims));
nStims   = numel(bins);

for f = fieldnames(data.spikecounts)'; f=f{1};
   n.(f) = size(data.spikecounts.(f),1);
   
   if strcmpi(data.protocol,'trans92')
      tun.(f).bins  = [reshape_sloppy(bins,12)'; NaN(7,8)];
      tun.(f).bins(:,5:8) = [[49:51 53:60 62:65 67:70];
                             [73:76 78:85 87:89 91:92 NaN(1,2)];
                             [52 61 66 71 NaN(1,15)];
                             [72 77 86 90 NaN(1,15)]]';
   else
      tun.(f).bins  = bins;
   end
   tun.(f).respmean = NaN(size(tun.(f).bins));
   tun.(f).respstd  = NaN(size(tun.(f).bins));
   tun.(f).protocol = data.protocol;

   tun.(f) = repmat(tun.(f),[n.(f) 1]);
   for j = 1:n.(f)
      for i = 1:nStims
         ix = data.stims==bins(i);
         binix = tun.(f)(j).bins==bins(i);
         tun.(f)(j).respmean(binix) = mean(data.spikecounts.(f)(j,ix ...
                                                        & data.valid),2);
         tun.(f)(j).respstd( binix) =  std(data.spikecounts.(f)(j,ix ...
                                                        & data.valid),[],2);
      end
      %% TODO: populate this
      tun.(f)(j).spontmean   = NaN;
      tun.(f)(j).spontstd    = NaN;
   end
end


return

