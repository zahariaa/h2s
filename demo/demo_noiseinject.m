% demo_variability

seed = rng(0);
SIMPLICES = false;%true;
another3  = false;
CONSTRAINT= false;
FIXRADII  = true;
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
   
      [model,high] = hypersphere2sphere(points,categories,[],dimLow);
      model = SetOfHyps(model,groundtruth);
      model.sig = model.significance(points,1000);
      %keyboard
      testhi = SetOfHyps('estimate',points,categories);%,100).merge;
      [~,~,testhi.error,testhi.msflips] = testhi.stress(groundtruth);
      testhi.sig = testhi.significance(points);
   
      %axtivate(fh.a.h(i*2*nconditions+j));  model.show;
      for MDS_INIT = [false true]
         testlo = testhi.h2s([dimLow 0 ninits],[FIXRADII CONSTRAINT MDS_INIT]);
         [~,~,testlo.error,testlo.msflips] = testlo.stress(groundtruth);
         testlo.sig = testlo.significance(points);
      
         axtivate(fh.a.h((i*2+double(MDS_INIT))*nconditions+j)); testlo.show;
         drawnow;
      end
   end
   %[groundtruth.overlap;model.overlap;testlo.overlap]
   i = i+1;
end

