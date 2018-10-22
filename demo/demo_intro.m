% demo_intro: generates introductory explanation figure

n     = 40;
nCats = 4;
dotsz = 6;

%% Simulate data
points = randnball(nCats*n,3);
% % line of spheres
% for i = 2:4
%    points((i-1)*n+1:i*n,1) = points((i-1)*n+1:i*n,1) + i-1;
% end
% tetrahedron of spheres
v = [1 0 -1/sqrt(2); -1 0 -1/sqrt(2); 0 1 1/sqrt(2); 0 -1 1/sqrt(2)];
for i = 1:4
   points((i-1)*n+1:i*n,:) = points((i-1)*n+1:i*n,:) + repmat(v(i,:),[n 1]);
end

categories = Categories(n*ones(1,4));

%% PLOT
planelim = 2;
for itype = 1:-1:0
   nCats = size(points,1)/n - itype;
   fh = newfigure(sprintf('intro%u',itype),[1 3]);
   axtivate(1)
   for i = 1:nCats
      meanPoint = mean(points((i-1)*n+1:i*n,:));
      plot3(points((i-1)*n+1:i*n,1),points((i-1)*n+1:i*n,2),...
            points((i-1)*n+1:i*n,3),'wo','MarkerSize',dotsz,...
            'MarkerFaceColor',categories.colors(i,:));
      for j = 1:n
         plot3([meanPoint(1) points((i-1)*n+j,1)],...
               [meanPoint(2) points((i-1)*n+j,2)],...
               [meanPoint(3) points((i-1)*n+j,3)],'-','Color',categories.colors(i,:))
      end
      plot3(meanPoint(1),meanPoint(2),meanPoint(3),...
            'wo','MarkerSize',dotsz*1.5,'MarkerFaceColor',categories.colors(i,:))
   end
   draw3Daxes([0 0 0],[-1 1 -1 1 -1 1]*planelim); axis vis3d off; view(40,22); rotate3d;

   axtivate(2) 
   sh = draw3dEllipsoid(v(1:nCats,:),[],categories.colors(1:nCats,:),[],1/3);
   draw3Daxes([0 0 0],[-1 1 -1 1 -1 1]*planelim); view(40,22);
   match3D(fh.a.h(1),fh.a.h(2));
   set(fh.f,'Renderer','openGL');
   axis(fh.a.h(3),'off');
   set(sh,'AmbientStrength',0,'SpecularStrength',1,'DiffuseStrength',1);
   camlight('left'); lighting phong

%% PLOT PLANE
   if itype==1
      % Compute intersection plane coefficients
      m = rref([v(1:3,:) ones(3,1)]);
      m = m(:,end);
      m = -m(1:2)/m(3);
      b = mean(v(1:3,:));
      % Draw plane
      [X,Y] = meshgrid([-planelim planelim],[-planelim planelim]);
      X = X([1 2 4 3 1]);
      Y = Y([1 2 4 3 1]);
      % limit plane extent to within axis limits (so plane outline is visible)
      X(X>0) = min(planelim,(planelim-b(3))/m(1)); X(X<0) = max(-planelim,(-planelim-b(3))/m(1));
      Y(Y>0) = min(planelim,(planelim-b(3))/m(2)); Y(Y<0) = max(-planelim,(-planelim-b(3))/m(2));
      plane = patch(X+b(1),Y+b(2),m(1)*X+m(2)*Y+b(3),zeros(size(X)),...
                    'EdgeColor','k','LineWidth',1,'FaceColor',[0.5 0.5 0],'FaceAlpha',0.5)
      % scaling for intersection circles
      mx = sqrt(1/(1+m(1)^2));
      my = sqrt(1/(1+m(2)^2));
      % Plot intersection circles
      t = linspace(0,2*pi,100);
      for i = 1:nCats
         circ(i) = plot3(mx*cos(t)+v(i,1),my*sin(t)+v(i,2),v(i,3)+mx*m(1)*cos(t)+my*m(2)*sin(t),...
                         'Color',categories.colors(i,:),'LineWidth',2);
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
      plot3(v(i,1),v(i,2),v(i,3),'wo','MarkerSize',dotsz,...
            'MarkerFaceColor',categories.colors(i,:));
      plot3(v(i,1)*[1 1]+[0 1],v(i,2)*[1 1],v(i,3)*[1 1],...
            'Color',categories.colors(i,:))
   end
   
   axtivate(3)
   model = hypersphere2sphere(points(1:(n*nCats),:),categories.select(1:nCats),[],2);
   showCirclesModel(model,[],[]);
   for i = 1:nCats
      plot(model.centers(i,1),model.centers(i,2),'wo','MarkerSize',dotsz,...
           'MarkerFaceColor',categories.colors(i,:))
      plot(model.centers(i,1)*[1 1]+[0 1],model.centers(i,2)*[1 1],...
           'Color',categories.colors(i,:))
   end
   draw3Daxes;
   papsz = [5.75 5.75/3];
   set(fh.f,'Renderer','painters','PaperUnits','inches','PaperSize',papsz,...
       'PaperPosition',[0.01*papsz papsz],'PaperPositionMode','manual');
   printFig;
end

%% TEST NEW STRESS FUNCTION
orig = SetOfHyps(v',ones(1,4));
low  = SetOfHyps(model);
low.error = low.stress(orig);
new  = orig.h2s
% Plot
fh = newfigure([1 2],'compare'); low.show(fh.a.h(1));
                                 new.show(fh.a.h(2));
set(fh.f,'Renderer','openGL')

