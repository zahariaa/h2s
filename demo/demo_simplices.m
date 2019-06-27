% demo_simplices


ds = [3:5 3];
another3 = false; i = 0;

fh = newfigure('n-simplices',[numel(ds) 2]);

for d = ds
   i = i+1;
   dimLow = 2;%min(d-1,3);
   if     d==3 &&  another3
      groundtruth = SetOfHyps([nsimplex(d)';nsimplex(d)'+3],ones(2*d+2,1));
   else
      another3    = true;
      groundtruth = SetOfHyps(nsimplex(d)',ones(d+1,1));
   end
   [points,categories] = randnsimplex_of_nballs(groundtruth.centers,50);
   [model,high] = hypersphere2sphere(points,categories,[],dimLow);
   model = SetOfHyps(model,groundtruth);
   %keyboard
   testhi = SetOfHyps('estimate',points,categories);
   testhi.error = testhi.stress(groundtruth);
   testlo = testhi.h2s(dimLow);
   testlo.error = testlo.stress(groundtruth);

   axtivate(fh.a.h((i-1)*2+1));  model.show;
   axtivate(fh.a.h((i-1)*2+2)); testlo.show;
   drawnow;
   [groundtruth.overlap;model.overlap;testlo.overlap]
end

