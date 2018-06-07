function TEST_hypersphere2sphere
% TEST_hypersphere2sphere

% hs2s tests: comparison to MDS
%       - two touching equal-radius hyperspheres in 1-10 dimensions
%       - two concentric hyperspheres of different radius in 1-10 dimensions
%       - larger hypersphere enclosing smaller one, touching in one surface point
%       - two intersecting hypershperes
%       - 3 spheres (two intersecting) -> 2D circles
%       - 5 hyperspheres in 1-100 dimensions
%       - as previous, but illustrate effect of sample size changes



%% general control variables
scenarios = [1 3:6]; %[1 2 3];
dimLow    = 2; % 2=showCircles, 3=showSpheres
nScen     = numel(scenarios);


PLOTALL = false;
psFilespec = '';
pdfFilespec = '';


%% initialize combo figure (scenarios 1, 3-6, for 3D and 200D); h2s in 2D & 3D
fh = newfigure([nScen 9],'whatsthescenario');
iCombo = 0;

nsDim = [1 2 3 5 10 20 40 200];

scenarios = {'scenario 1: two touching equal-radius hyperspheres in 1-200 dimensions';
             'scenario 2: two identical hyperspheres in 1-200 dimensions';
             'scenario 3: two concentric hyperspheres of different radius in 1-200 dimensions';
             'scenario 4: larger hypersphere enclosing smaller one, touching in one surface point';
             'scenario 5: two intersecting hyperspheres';
             'scenario 6: two intersecting hyperspheres with different numbers of points';
             'scenario 7: 2 overlap hyperellipsoids w rand covar 1-200 dims';
             'scenario 8: 2 adj hyperellipsoids w same rand covar 1-200 dims';
             'scenario 9: 2 overlap hyperellipsoids w rand covar 1-200 dims'};
centers = {[2 0],[0 0],[0 0],[0 1],[1 0],[1 0],[0 0],[1 -1],[1 -1]};
radii   = {[1 1],[1 1],[2 1],[2 1],[1 1],[1 1],[1 1],[1 1],[1 1]};
nPtsPerCat = {[40 40],[100 100],[50 50],[100 100],[50 50],[20*5 20],[40 40],[40 40],[50 50*5]};

%% MAIN LOOP
for s = scenarios
   if PLOTALL, fig(s)=newfigure(scenarios{s}); end
   for nDimI = 1: numel(nsDim)
      nDim = nsDim(nDimI);
      if s < 7 % UNIFORM
         points = [randsphere(nPtsPerCat{s}(1),nDim,radii{s}(1))
                   randsphere(nPtsPerCat{s}(2),nDim,radii{s}(2))];
         points(1:nPtsPerCat{s}(1)    ,1) = points(1:nPtsPerCat{s}(1)    ,1) + centers{s}(1);
         points(nPtsPerCat{s}(1)+1:end,1) = points(nPtsPerCat{s}(1)+1:end,1) + centers{s}(2);
      else    % GAUSSIAN
         SIG  = rand(nDim,nDim,2);
         for i = 1:2;   SIG(:,:,i) = SIG(:,:,i)*SIG(:,:,i)';   end
         if s>=8, SIG(:,:,2) = SIG(:,:,1); end
         points = nan(sum(nPtsPerCat{s}),nDim);
         points(1:nPtsPerCat{s}(1),:)     = mvnrnd(centers{s}(1)*ones(nDim,1),SIG(:,:,1),nPtsPerCat{s}(1));
         points(nPtsPerCat{s}(1)+1:end,:) = mvnrnd(centers{s}(2)*ones(nDim,1),SIG(:,:,2),nPtsPerCat{s}(1));
      end
      categories = Categories(nPtsPerCat{s});
      if PLOTALL
      figPanelSpec = [fig(s) 4 8 1+(nDimI-1)*4];
      titleStr = any2str(nPtsPerCat{s}(1), ' points/cat. in ',nDim,' dim.');
      HS2SandMDS(points,categories,figPanelSpec,titleStr,dimLow)
      end
      iCombo = updateComboPlot(fh,iCombo,nDim,points,categories,dimLow,centers{s},radii{s});
   end
end

%% Finish up combo figure
%set(fh.a.h(1:9:end),'CameraPosition',[0 0 -33]); % make 3D views consistent
papsz = [6.5 4];
set(fh.f,'PaperUnits','inches','PaperSize',papsz,'PaperPosition',[0.01*papsz papsz],'PaperPositionMode','manual');
set(findobj(fh.f,'type','line'),'MarkerSize',3); % make dots slightly smaller
subplotResize([],0.01,0.01);
set(fh.f,'Renderer','painters');   printFig;
delete(fh.a.h(2:9:end));
set(fh.f,'Renderer','openGL');     printFig(fh.f,[],'png',1200);

%% Add 3D & 2D renders to "Combo plot" for paper
function iCombo = updateComboPlot(fh,iCombo,nDim,points,categories,dimLow,centers,radii)
if     nDim==3
   iCombo = iCombo+1;
   % format abbreviated radii and centers into full 3D vectors/matrices
   locations  = zeros(numel(centers),3); locations(:,1) = centers(:);
   covariance = arrayfun(@(r) r*eye(3),radii.^2,'UniformOutput',false)';

   axtivate(fh.a.h(1+(iCombo-1)*9));
   draw3dEllipsoid(locations,covariance,categories.colors,[],1/3);
   view([180 90]);

   ax = fh.a.h((2:5)+(iCombo-1)*9);
   HS2SandMDS(points,categories,ax,[],dimLow);
elseif nDim==200
   ax = fh.a.h((6:9)+(iCombo-1)*9);
   HS2SandMDS(points,categories,ax,[],dimLow);
end

