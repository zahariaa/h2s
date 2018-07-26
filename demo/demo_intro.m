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
fh = newfigure('intro',[1 3]);
axtivate(1)
for i = 1:nCats
   plot3(points((i-1)*n+1:i*n,1),points((i-1)*n+1:i*n,2),...
         points((i-1)*n+1:i*n,3),'wo','MarkerSize',dotsz,...
         'MarkerFaceColor',categories.colors(i,:));
end
draw3Daxes([0 0 0],[-1 1 -1 1 -1 1]*2); axis vis3d off; view(40,22); rotate3d;

axtivate(2) 
sh = draw3dEllipsoid(v,[],categories.colors,[],1/3);
draw3Daxes([0 0 0],[-1 1 -1 1 -1 1]*2); view(40,22);
match3D(fh.a.h(1),fh.a.h(2));
set(fh.f,'Renderer','openGL');
axis(fh.a.h(3),'off');
set(sh,'AmbientStrength',0,'SpecularStrength',1,'DiffuseStrength',1);
camlight('left'); lighting phong
% export spheres in png
subplotResize(fh.a.h,[],0.01); printFig([],[],'png',200);
delete(sh); % delete spheres for pdf

for i=1:nCats % Plot centers and radii
   plot3(v(i,1),v(i,2),v(i,3),'wo','MarkerSize',dotsz,...
         'MarkerFaceColor',categories.colors(i,:));
   plot3(v(i,1)*[1 1]+[0 1],v(i,2)*[1 1],v(i,3)*[1 1],...
         'Color',categories.colors(i,:))
end

axtivate(3)
model = hypersphere2sphere(points,categories,[],2);
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


%% TEST NEW STRESS FUNCTION
orig = Hyperspheres(v',ones(1,4));
low  = Hyperspheres(model);
low.error = low.stress(orig);
new  = orig.h2s
% Plot
fh = newfigure([1 2],'compare'); low.show(fh.a.h(1));
                                 new.show(fh.a.h(2));
set(fh.f,'Renderer','openGL')

