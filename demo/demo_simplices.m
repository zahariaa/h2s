% demo_simplices


ds = 3:5;

fh = newfigure('n-simplices',[numel(ds) 2]);

for d = ds
   dimLow = 2;%min(d-1,3);
   groundtruth = Hyperspheres(nsimplex(d)',ones(d+1,1));
   [points,categories] = randnsimplex_of_nballs(d,50);
   [model,high] = hypersphere2sphere(points,categories,[],dimLow);
   model = Hyperspheres(model,groundtruth);
   testhi = Hyperspheres('estimate',points,categories);
   testhi.error = testhi.stress(groundtruth);
   testlo = testhi.h2s(dimLow);

   axtivate((find(ds==d)-1)*2+1); model.show(dimLow);
   axtivate((find(ds==d)-1)*2+2); testlo.show(dimLow);
   drawnow;
end

