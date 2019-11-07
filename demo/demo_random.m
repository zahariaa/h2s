% demo_random

seed = rng(0);
SIMPLICES = false; %true;
another3  = false;

if SIMPLICES
   ds = [3 5 7 3];
else
   ds = [3 4 6 9];
end

nconditions = 5;
fh = newfigure('n-random',[numel(ds) nconditions]);
set(fh.f,'Units','normalized');

i = 0;
for d = ds
   i = i+1;
   j = 1;
   dimLow = min(d-1,3);
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
   [model,high] = hypersphere2sphere(points,categories,[],dimLow);
   model = SetOfHyps(model,groundtruth);
   model.sig = model.significance(points,1000);
   %keyboard
   testhi = SetOfHyps('estimate',points,categories,groundtruth);%,100).merge;
   [testhi.sig,sechi] = testhi.significance(points);
   axtivate(fh.a.h((i-1)*nconditions+j));  model.show;

   for MDS_INIT = true%[false true]
      for FIXRADII = false%[false true]
         j = j+1;
         testlo = testhi.h2s(dimLow,[FIXRADII MDS_INIT],groundtruth);
         [testlo.sig,seclo] = testlo.significance(points,1000);
      
         testlo.show(      fh.a.h((i-1)*nconditions+j  ));
         testlo.showValues(fh.a.h((i-1)*nconditions+j+1),testhi);
         testhi.showSig(   fh.a.h((i-1)*nconditions+j+2),'legend');
         testhi.showSig(   fh.a.h((i-1)*nconditions+j+3),sechi);
         drawnow;
      end
   end
   %[groundtruth.overlap;model.overlap;testlo.overlap]
end

