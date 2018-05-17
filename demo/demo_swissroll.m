function demo_swissroll

dimLow = 2;

%% set up swiss roll
rng(10); % for reproducibility
N = 150;
noise = 0.1;
t = 3*pi/2 * (1 + 2*rand(N,1));
r = 5 * rand(N,1);

%% set up points and categories
ca.labels  = {'category 1','category 2'};
ca.colors  = [0.8 0 0; 0 0 0];
ca.vectors = logical(blkdiag(true(N,1),true(N,1)));

swissroll  = { @(t,r) [t.*cos(t   ), r, t.*sin(t   )];
               @(t,r) [t.*cos(t+pi), r, t.*sin(t+pi)] };

points = noise*randn(2*N,3);
points(ca.vectors(:,1),:) = points(ca.vectors(:,1),:) + swissroll{1}(t,r);
points(ca.vectors(:,2),:) = points(ca.vectors(:,2),:) + swissroll{2}(t,r);

%% PLOT
figure;hold on;
plot3(points(ca.vectors(:,1),1),points(ca.vectors(:,1),2),points(ca.vectors(:,1),3),'ko')
plot3(points(ca.vectors(:,2),1),points(ca.vectors(:,2),2),points(ca.vectors(:,2),3),'ro')
axis vis3d off;view(-6,2);draw3Daxes

%% H2S
h=figure(2010);clf;
set(h,'Name','H2S swiss roll');
i=1;
figPanelSpec = [h 2 4 1+(i-1)*4];
titleStr = sprintf('%u points/cat',N);
HS2SandMDS(points,ca,figPanelSpec,titleStr,dimLow)


%% siiiiiiine waaaaaves
ca.labels  = {'category 1','category 2','category 3','category 4'};
ca.colors  = [0 0 0; 0.8 0 0; 0 0.8 0; 0 0 0.8];
ca.vectors = logical(blkdiag(true(N,1),true(N,1),true(N,1),true(N,1)));

sinewaves  = { @(t,r) [r, t, 3*cos(t)          ];
               @(t,r) [r, t, 3*cos(t)+5        ]; 
               @(t,r) [r, t, 3*cos(t)+2*r+2*t  ]; 
               @(t,r) [r, t, 3*cos(t)+2*r+2*t+5] }; 

points = noise*randn(4*N,3);
for i = 1:numel(ca.labels)
   points(ca.vectors(:,i),:) = points(ca.vectors(:,i),:) + sinewaves{i}(t,r);
end

%% PLOT
figure;hold on;
set(gcf,'Name','sine waves');
for i = 1:numel(ca.labels)
   plot3(points(ca.vectors(:,i),1),points(ca.vectors(:,i),2),...
         points(ca.vectors(:,i),3),'o','MarkerFaceColor',ca.colors(i,:),...
         'MarkerEdgeColor','w')
end
axis vis3d off;view(88,6);draw3Daxes

%% H2S
i=2;
figPanelSpec = [h 2 4 1+(i-1)*4];
titleStr = sprintf('%u points/cat',N);
HS2SandMDS(points,ca,figPanelSpec,titleStr,dimLow)


