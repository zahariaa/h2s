% QT_uniformWithinHypersphereSummaryStats

%% control variables
nsDim=1:4;
nsDim=[1 2 3 6 15 40 100];
nsDim=[1 2 3 6 10 20 40 500 1000];

lw=2;




%% simulated densities
n=1e5;
i=1;
cols=colorScale([.85 .85 .85; 0 0 0],numel(nsDim));
stds=[];
legendStrings=[];

h=figure(10); clf; set(h,'Color','w');
subplot(2,2,1);

for nDim=nsDim
    %points=randn(n,nDim);
    points=randsphere(n,nDim,1);

    radii = sqrt(sum(points.^2,2));
    stds(i) = std(radii);
    mean_dist2center(i) = mean(radii);
    
    subplot(2,2,1);
    [counts,xout]=hist(radii,50);
    plot(xout,counts,'Color',cols(i,:),'LineWidth',lw); hold on;
    
    points_centered = points - repmat(mean(points,1),[n 1]);
    radii_centered = sqrt(sum(points_centered.^2,2));
    
    subplot(2,2,2);
    [counts,xout]=hist(radii_centered,50);
    plot(xout,counts,'Color',cols(i,:),'LineWidth',lw); hold on;
    
    n2=400;
    p2 = points(1:n2,:);
    dists = pdist(p2,'Euclidean');
    mean_dist(i) = mean(dists);
    
    subplot(2,2,3);
    [counts,xout]=hist(dists,50);
    plot(xout,counts,'Color',cols(i,:),'LineWidth',lw); hold on;
    
    legendStrings{i}=[num2str(nDim),' dim.'];
    i=i+1;
end
subplot(2,2,1); axis([0 1 0 n/5]); axis square;
xlabel('radius'); ylabel('count');
title('\bfradii from true center');

subplot(2,2,2); axis([0 1 0 n/5]); axis square;
xlabel('radius'); ylabel('count');
title('\bfradii from center estimated as mean');
legend(legendStrings,'Location','NorthWest');

subplot(2,2,3); axis([0 2 0 n2^2/30]); axis square;
xlabel('distance between two points'); ylabel('count');
title('\bfdistances between points');

subplot(2,2,4);
plot(nsDim,stds,'k');
xlabel('number of dimensions');
ylabel('standard deviation of radii');



%% plot expected distance to center and distance between two points
% as a function of dimensionality
h=figure(20); clf; set(h,'Color','w');
lw=3;
ms=20;

dom = 0:max(nsDim);

subplot(2,1,1); hold on;
%plot(dom,sqrt(dom),'Color',[0.7 0.7 0.7],'LineWidth',lw);
plot(nsDim,mean_dist2center,'k.','MarkerSize',ms);
xlabel('number of dimensions');
ylabel('expected distance to center');

subplot(2,1,2); hold on;
%plot(dom,sqrt(dom*2),'Color',[0.7 0.7 0.7],'LineWidth',lw);
plot(nsDim,mean_dist,'k.','MarkerSize',ms);
xlabel('number of dimensions');
ylabel('expected distance between two points');



    