% demo_books: h2s on book word frequency

datafile = 'books';
load(['~' filesep datafile '.mat']);

for i = 1:numel(titles)
   titles{i}(strfind(titles{i},' ')) = [];
end

chapix  = [0 cumsum(chapix)];
vectors = false(size(freq,1),numel(titles));
newfreq = [];
for i = 1:numel(titles)
   newfreq = blkdiag(newfreq,squeeze(freq(1+chapix(i):chapix(i+1),:,:)));
   vectors(1+chapix(i):chapix(i+1),i) = true;
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

cats = Categories(vectors,titles);

hi = SetOfHyps('estimate',newfreq,cats.select([1:2 5:7]),1);

figure;hi.h2s(3).show;cats.select([1:2 5:7]).legend;


