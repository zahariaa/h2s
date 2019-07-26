% demo_random

seed = rng(0);
SIMPLICES = true;
another3  = false;

ds = [3 3 4 4 6 6 9 9];%[3:5 3];

nconditions = 5;
fh = newfigure('n-random',[numel(ds) nconditions]);

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
   %keyboard
   testhi = SetOfHyps('estimate',points,categories);%,100).merge;
   testhi.error = testhi.stress(groundtruth);
   axtivate(fh.a.h((i-1)*nconditions+j));  model.show;

   for CONSTRAINT = [false true]
      for FIXRADII = [false true]
         j = j+1;
         testlo = testhi.h2s(dimLow,[FIXRADII CONSTRAINT]);
         testlo.error = testlo.stress(groundtruth);
      
         axtivate(fh.a.h((i-1)*nconditions+j)); testlo.show;
         drawnow;
      end
   end
   %[groundtruth.overlap;model.overlap;testlo.overlap]
end

