% demo_noiseinject

seed = rng(1);
SIMPLICES = true;
another3  = false;
CONSTRAINT= false;
FIXRADII  = true;
ninits    = 5;

if SIMPLICES
   ds = [3 5 7 3];
else
   ds = [3 4 6 9];
end

% fh = newfigure('n-random',[numel(ds) nconditions]);
% set(fh.f,'Units','normalized');

for d = ds
   if SIMPLICES
      if     d==3 &&  another3
         groundtruth = SetOfHyps([nsimplex(d)';nsimplex(d)'+3],ones(2*d+2,1));
      else
         another3    = true;
         groundtruth = SetOfHyps(nsimplex(d)',ones(d+1,1));
      end
   else
      groundtruth = SetOfHyps(2*rand(d+1,d),rand(d+1,1));
   end
   [points,categories] = randnsimplex_of_nballs(groundtruth.centers,50,...
                                                groundtruth.radii);

   hi = SetOfHyps('estimate',points,categories);%,100).merge;
   rfixed = hi.h2s(true); % MDS_INIT & FIXEDRADII=true
   rfixed.sig = rfixed.significance(points);

   rfree  = rfixed.h2s(hi,false);
   rfree.sig = rfree.significance(points);

   % All error
   [sum(rfixed.error) sum(rfree.error)]
   % center errors only
   nc2 = nchoosek(numel(hi.radii),2);
   [sum(rfixed.error(1:nc2)) sum(rfree.error(1:nc2))]

   %axtivate(fh.a.h((i*2+double(MDS_INIT))*nconditions+j)); lo.show;
   %drawnow;
end

