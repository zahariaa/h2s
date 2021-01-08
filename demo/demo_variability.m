% demo_variability

seed = rng(1);
SIMPLICES = false;%true;
another3  = false;
ninits    = 5;

if SIMPLICES
   ds = [3 5 7 3];
else
   ds = [3 4 6 9];
end

nconditions = 5;
fh = newfigure('n-random',[numel(ds)*2 nconditions]);
set(fh.f,'Units','normalized');

i = 0;
for d = ds
   dimLow = min(d-1,3);

   for j = 1:nconditions
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
   
      %keyboard
      testhi = SetOfHyps.estimate(points,categories,groundtruth,1000).meanAndMerge;
   
      for MDS_INIT = {'randinit' 'mdsinit'}
         testlo = testhi.h2s([dimLow ninits],MDS_INIT{1});
         testlo = testlo.stressUpdate(groundtruth);
         testlo = testlo.significance(points,1000);
      
         axtivate(fh.a.h((i*2+double(strcmpi(MDS_INIT,'mdsinit')))*nconditions+j));
         testlo.show;
         drawnow;
      end
   end
   i = i+1;
end

