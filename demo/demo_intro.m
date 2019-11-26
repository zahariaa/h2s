% demo_intro: generates introductory explanation figure

n     = 50;
nCats = 4;
dotsz = 6;

%% Simulate data
% % line of spheres
% for i = 2:4
%    points{1}((i-1)*n+1:i*n,1) = points{1}((i-1)*n+1:i*n,1) + i-1;
% end
% tetrahedron of spheres
v{1} = [9 0 0; [9 9]/sqrt(2) 0; [4 4]/sqrt(2) 0; 0 0 0] - [1;1;1;1]*[6 0 0];
r = {1:4,ones(1,6)/sqrt(3/2)};

%%TODO: generalize simplex scale factor to have 0 overlap
%d=3;hi = SetOfHyps(nsimplex(d)*sqrt(3/2),ones(d+1,1));hi.margins
v{2} = nsimplex(3)';
v{2} = v{2}([1:3 2:4],:) + [[3;3;3;-3;-3;-3] zeros(6,2)];

points = cellfun(@(x) randnball(n*numel(x),3),r,'UniformOutput',false);
for j = 1:numel(r)
   for i = 1:numel(r{j})
      points{j}((i-1)*n+1:i*n,:) = r{j}(i)*points{j}((i-1)*n+1:i*n,:) ...
                                   + repmat(v{j}(i,:),[n 1]);
   end
end


%% PLOT
planelim = 2;
for itype = 0%1%1:-1:0
   orig = SetOfHyps(v{itype+1},r{itype+1});
   orig.categories = Categories(n*ones(1,numel(r{itype+1})));

   nCats = size(points{itype+1},1)/n;
   fh = newfigure(sprintf('intro%u',itype),[1 3]);
   axtivate(1)
   for i = 1:nCats
      plot3(points{itype+1}((i-1)*n+1:i*n,1),points{itype+1}((i-1)*n+1:i*n,2),...
            points{itype+1}((i-1)*n+1:i*n,3),'wo','MarkerSize',dotsz,...
            'MarkerFaceColor',orig.categories.colors(i,:));
   end
   draw3Daxes([0 0 0],5*[-1 1 -1 1 -1 1]*planelim); axis equal vis3d off; view(40,22); rotate3d;

   axtivate(2) 
   sh = draw3dEllipsoid(v{itype+1},...
            arrayfun(@(x) eye(3)*x^2,r{itype+1},'UniformOutput',false),...
            orig.categories.colors(1:nCats,:),[],0.4);
   draw3Daxes([0 0 0],5*[-1 1 -1 1 -1 1]*planelim); view(40,22);
   match3D(fh.a.h(1),fh.a.h(2));
   set(fh.f,'Renderer','openGL');
   axis(fh.a.h(3),'off');
   set(sh,'AmbientStrength',0,'SpecularStrength',1,'DiffuseStrength',1);
   camlight('left'); lighting phong

%% PLOT PLANE
   if itype==0
      % Compute intersection plane coefficients
      m = rref([v{itype+1}(1:3,:) ones(3,1)]);
      m = m(:,end);
      m = -m(1:2)/m(3);
      b = mean(v{itype+1}(1:3,:));
      % Draw plane
      [X,Y] = meshgrid([-planelim planelim],[-planelim planelim]);
      X = X([1 2 4 3 1]);
      Y = Y([1 2 4 3 1]);
      % limit plane extent to within axis limits (so plane outline is visible)
      X(X>0) = min(planelim,(planelim-b(3))/m(1)); X(X<0) = max(-planelim,(-planelim-b(3))/m(1));
      Y(Y>0) = min(planelim,(planelim-b(3))/m(2)); Y(Y<0) = max(-planelim,(-planelim-b(3))/m(2));
      plane = patch(X+b(1),Y+b(2),m(1)*X+m(2)*Y+b(3),zeros(size(X)),...
                    'EdgeColor','k','LineWidth',1,'FaceColor',[0.5 0.5 0.5],'FaceAlpha',0.5);
      % scaling for intersection circles
      mx = sqrt(1/(1+m(1)^2));
      my = sqrt(1/(1+m(2)^2));
      % Plot intersection circles
      t = linspace(0,2*pi,100);
      for i = 1:nCats
         circ(i) = plot3(v{itype+1}(i,1) + r{itype+1}(i)*mx*cos(t),...
                         v{itype+1}(i,2) + r{itype+1}(i)*my*sin(t),...
                         v{itype+1}(i,3) + r{itype+1}(i)*mx*cos(t)*m(1) ...
                                         + r{itype+1}(i)*my*sin(t)*m(2),...
                         'Color',orig.categories.colors(i,:),'LineWidth',2);
      end
   else circ = []; plane = [];
   end
   % %DEBUG
   % axis([-4 4 -4 4 -4 4])
   % keyboard
   % export spheres in png
   subplotResize(fh.a.h,[],0.01); printFig([],[],'png',200);
   delete([sh(:);circ(:);plane]);   % delete spheres, intersection circles, and plane for pdf

   for i=1:nCats % Plot centers and radii
      plot3(v{itype+1}(i,1),v{itype+1}(i,2),v{itype+1}(i,3),'wo','MarkerSize',dotsz,...
            'MarkerFaceColor',orig.categories.colors(i,:));
      plot3(v{itype+1}(i,1)*[1 1]+[0 r{itype+1}(i)],...
            v{itype+1}(i,2)*[1 1],v{itype+1}(i,3)*[1 1],...
            'Color',orig.categories.colors(i,:))
   end
   % Plot lines connecting centers
   for i = 1:nCats
      for j = i+1:nCats
         plot3(v{itype+1}([i j],1),v{itype+1}([i j],2),v{itype+1}([i j],3),'-','Color',[0.5 0.5 0.5])
      end
   end
   
   axtivate(3)
   model = hypersphere2sphere(points{itype+1}(1:(n*nCats),:),orig.categories.select(1:nCats),[],2);
   [~,~,model.error] = SetOfHyps(model).stress(orig);
   showCirclesModel(model,[],[]);
   for i = 1:nCats
      plot(model.centers(i,1),model.centers(i,2),'wo','MarkerSize',dotsz,...
           'MarkerFaceColor',orig.categories.colors(i,:))
      plot(model.centers(i,1)*[1 1]+[0 model.radii(i)],...
           model.centers(i,2)*[1 1],'Color',orig.categories.colors(i,:))
   end
   % Plot lines connecting centers
   for i = 1:nCats
      for j = i+1:nCats
         plot(model.centers([i j],1),model.centers([i j],2),'-','Color',[0.5 0.5 0.5])
      end
   end
   draw3Daxes;
   papsz = [5.75 5.75/3];
   set(fh.f,'Renderer','painters','PaperUnits','inches','PaperSize',papsz,...
       'PaperPosition',[0.01*papsz papsz],'PaperPositionMode','manual');
   printFig;
end

%% TEST NEW STRESS FUNCTION
low  = SetOfHyps(model);
[~,~,low.error] = low.stress(orig);
new  = low.h2s(orig,2)%orig.h2s
% Calculate significance
[orig.sig,sec] = orig.significance(points{itype+1},10000);
% Plot
low.sig = orig.sig;
new.sig = orig.sig;
fh = newfigure([3 2],'compare'); low.show(fh.a.h(1));
                                 new.show(fh.a.h(2));
set(fh.f,'Renderer','openGL')
orig.showSig(fh.a.h(5),'legend');
orig.showSig(sec,fh.a.h(6));
matchy(fh.a.h(5:6))

low.showValues(fh.a.h(3))
new.showValues(low,fh.a.h(4))

% Print
set(fh.f,'Renderer','painters')
printFig(fh.f,[],'eps')

fh = newfigure([1 4],'stats');
new.show(fh.a.h(1));
new.showValues(low,fh.a.h(2))
orig.showSig(fh.a.h(3));
orig.showSig(sec,fh.a.h(4));
matchy(fh.a.h(2:4))
% Print
set(fh.f,'Renderer','painters')
printFig(fh.f,[],'eps')

