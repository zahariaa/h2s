% demo_random

%% USER-CHOSEN PARAMETERS
SEED      = rng(1);
SIMPLICES = true;
ALLOW3D   = false;
DBUGPLOT  = 'quiet'; %'verbose' %'debug' %'quiet'

if SIMPLICES
   ds = [3 5 7 3];
else
   ds = [3 4 6 9];
end

%% FIGURE SETUP
if SIMPLICES
   figtitle = sprintf('n-simplices_h2s%u',ALLOW3D+2);
else
   figtitle = sprintf('n-random_h2s%u'   ,ALLOW3D+2);
end
nconditions = 5;
another3  = false;
fh = newfigure(figtitle,[numel(ds) nconditions]);
set(fh.f,'Units','normalized');

%% FIGURE GENERATION LOOP
i = 0;
for d = ds
   i = i+1;
   j = 1;
   if ALLOW3D, dimLow = min(d-1,3);
   else        dimLow = 2;
   end
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
   model = Hypersphere.estimate(points,categories).h2s(dimLow);
   model = SetOfHyps(model).stressUpdate(groundtruth);
   % model = model.significance(points,1000);
   %keyboard
   testhi = SetOfHyps('estimate',points,categories,groundtruth,1000);
   axtivate(fh.a.h((i-1)*nconditions+j));  model.show;

   for MDS_INIT = {'mdsinit'} %{'mdsinit' 'randinit'}
      for alpha = {[1 1 1]} %{[1 1 1],[1 1 0]} %fix radii or not
         j = j+1;
         testlo = testhi.h2s(dimLow,MDS_INIT{1},DBUGPLOT,groundtruth,alpha);
         testlo = testlo.significance(points,1000);
      
         testlo.show(      fh.a.h((i-1)*nconditions+j  ));
         testlo.showValues(fh.a.h((i-1)*nconditions+j+1),testhi);
         testhi.showSig(   fh.a.h((i-1)*nconditions+j+[2 3]),'legend');
         drawnow;
      end
   end
   %[groundtruth.overlap;model.overlap;testlo.overlap]
end

%% PRINT
if ALLOW3D, printFig(fh,[],'png',450);
else        printFig(fh,[],'eps');
end

