function demo_swissroll

dimLow   = 2;
swlength = 5; % length of sine wave planes on axis orthogonal to wave

%% set up samples for swiss roll and sine waves
rng(10); % for reproducibility
N = 150;
noise = 0.1;
t = 3*pi/2 * (1 + 2*rand(N,1));
r = swlength * rand(N,1);

%% set up meshgrid for mesh versions of plots
[T,R] = meshgrid(3*pi/2 * linspace(1,3,125),linspace(0,swlength,swlength));

%% set up points and categories
ca.labels  = {'category 1','category 2'};
ca.colors  = [0 0 0; 0.8 0 0];
ca.vectors = logical(blkdiag(true(N,1),true(N,1)));

swissroll  = { @(t,r) [t.*cos(t   ), r, t.*sin(t   )];
               @(t,r) [t.*cos(t+pi), r, t.*sin(t+pi)] };

points = noise*randn(2*N,3);
points(ca.vectors(:,1),:) = points(ca.vectors(:,1),:) + swissroll{1}(t,r);
points(ca.vectors(:,2),:) = points(ca.vectors(:,2),:) + swissroll{2}(t,r);

%% PLOT
newfigure('swiss roll');
for i = 1:numel(ca.labels)
   srf = reshape(swissroll{i}(T(:),R(:)),[size(T) 3]);
   surf(srf(:,:,1),srf(:,:,2),srf(:,:,3),zeros(size(T(:))),...
        'FaceColor',ca.colors(i,:),'EdgeColor','none','FaceAlpha',0.5)
   plot3(points(ca.vectors(:,i),1),points(ca.vectors(:,i),2),...
         points(ca.vectors(:,i),3),'o','MarkerFaceColor',ca.colors(i,:),...
                                       'MarkerEdgeColor','w')
end
axis vis3d off;view(-6,2);draw3Daxes

%% H2S
fh = newfigure([2 4],'H2S swiss roll');
titleStr = sprintf('%u points/cat',N);
HS2SandMDS(points,ca,fh.a.h(1:4),titleStr,dimLow)


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
newfigure('sine waves');
for i = 1:numel(ca.labels)
    srf = reshape(sinewaves{i}(T(:),R(:)),[size(T) 3]);
    surf(srf(:,:,1),srf(:,:,2),srf(:,:,3),zeros(size(T(:))),...
        'FaceColor',ca.colors(i,:),'FaceAlpha',0.5)
   plot3(points(ca.vectors(:,i),1),points(ca.vectors(:,i),2),...
         points(ca.vectors(:,i),3),'o','MarkerFaceColor',ca.colors(i,:),...
                                       'MarkerEdgeColor','w')
end
axis vis3d off;view(88,6);draw3Daxes

%% H2S
titleStr = sprintf('%u points/cat',N);
HS2SandMDS(points,ca,fh.a.h(5:8),titleStr,dimLow)

%% MESH PLOTS
newfigure('meshes');
subplot(2,1,1);hold on;
for i = 1:2
   srf = reshape(swissroll{i}(T(:),R(:)),[size(T) 3]);
   surf(srf(:,:,1),srf(:,:,2),srf(:,:,3),zeros(size(T(:))),...
        'FaceColor',ca.colors(i,:),'EdgeColor','none','FaceAlpha',0.5)
end
axis vis3d off;view(-8,-6);draw3Daxes;
subplot(2,1,2);hold on;
for i = 1:numel(ca.labels)
    srf = reshape(sinewaves{i}(T(:),R(:)),[size(T) 3]);
    surf(srf(:,:,1),srf(:,:,2),srf(:,:,3),zeros(size(T(:))),...
        'FaceColor',ca.colors(i,:),'EdgeColor','none','FaceAlpha',0.5)
end
axis vis3d off;view(88,6);draw3Daxes;

