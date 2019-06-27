% demo_random

seed = rng(0);

ds = [3 3 4 4 5 5 8 8];%[3:5 3];

fh = newfigure('n-random',[numel(ds) 2]);

i=0;
for d = ds
   i = i+1;
   dimLow = min(d-1,3);
   groundtruth = SetOfHyps(2*rand(d+1,d),rand(d+1,1));

   [points,categories] = randnsimplex_of_nballs(groundtruth.centers,50,...
                                                groundtruth.radii);
   [model,high] = hypersphere2sphere(points,categories,[],dimLow);
   model = SetOfHyps(model,groundtruth);
   %keyboard
   testhi = SetOfHyps('estimate',points,categories);%,100).merge;
   testhi.error = testhi.stress(groundtruth);
   testlo = testhi.h2s(dimLow);
   testlo.error = testlo.stress(groundtruth);

   axtivate(fh.a.h((i-1)*2+1));  model.show;
   axtivate(fh.a.h((i-1)*2+2)); testlo.show;
   drawnow;
   %[groundtruth.overlap;model.overlap;testlo.overlap]
end

