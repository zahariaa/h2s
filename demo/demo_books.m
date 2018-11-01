% demo_books: h2s on book word frequency

datafile = 'books';
load(['~' filesep datafile '.mat']);

for i = 1:numel(titles)
   titles{i}(findstr(titles{i},' ')) = [];
end

chapix  = [1 cumsum(chapix)];
vectors = false(size(freq,1),numel(titles));
for i = 1:numel(titles)
   vectors(chapix(i):chapix(i+1),i) = true;
end

cats = Categories(vectors,titles);

hi = SetOfHyps('estimate',freq,cats.select([1:2 5:7]),1);

figure;hi.h2s(3).show;cats.select([1:2 5:7]).legend;

