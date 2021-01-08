% demo_books: h2s on book word frequency

%% Load
datafile = 'books';
load(['~' filesep datafile '.mat']);

freq = double(freq);
%% Strip spaces from legend text (why is this necessary?)
for i = 1:numel(titles)
   titles{i}(strfind(titles{i},' ')) = [];
end
%% Format vectors
vectors = false(size(freq,1),numel(titles));
for i = 1:numel(titles)
   vectors(1+chapix(i):chapix(i+1),i) = true;
end

%% Merge separate frequency matrices if given
if ndims(freq)==3
   newfreq = [];
   for i = 1:numel(titles)
      newfreq = blkdiag(newfreq,squeeze(freq(1+chapix(i):chapix(i+1),:,:)));
   end
   [words,wix] = sort(words);
   newfreq = newfreq(:,wix);
   
   for i = 1:numel(wix)-1
      if strcmp(words{i},words{i+1})
         newfreq(:,i  ) = sum(newfreq(:,i:i+1),2);
         newfreq(:,i+1) = NaN(chapix(end),1);
      end
   end
   words = words(~isnan(newfreq(1,:)));
   newfreq = newfreq(:,~isnan(newfreq(1,:)));
   freq = newfreq;
end
% %% Normalize if using CountVectorizer
% nw=sum(freq,2);
% freq(nw~=0,:) = freq(nw~=0,:)./repmat(nw(nw~=0),[1 numel(words)]);

%% COMPUTE
booksel = 1:6;%[1:3 6:7];
cats = Categories(vectors,titles);
hi   = SetOfHyps.estimate(freq,cats.select(booksel),1);
%% PLOT
figure; imagesc(freq); colorbar
figure; hi.h2s(3).show; cats.select(booksel).legend;
rotate3d


