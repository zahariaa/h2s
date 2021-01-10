% QT_GaussianDensityByRadiusInHiDim

%% control variables
nsDim=1:8;
nsDim=[1 2 3 6 15 40 100];
%nsDim=[1 2 3 6 10 20 40 500 1000];

lw=2;


%% analytical densities
r=0:0.01:200;
gaussian1D=1/(2*pi)*exp(-(r.^2/2));

cols=colourScale([.85 .85 .85; 0 0 0],numel(nsDim)); 
i=1;
h=figure(10); clf; set(h,'Color','w');
subplot(2,2,1);

for nDim=nsDim
    densityByR=gaussian1D.*(r.^(nDim-1));
    densityByR=densityByR./sum(densityByR)*sqrt(nDim);
    plot(r/sqrt(nDim),densityByR,'Color',cols(i,:),'LineWidth',lw); hold on;
    legendStrings{i}=[num2str(nDim),' dim.'];
    i=i+1;
end

axis([0 3 0 0.06]); axis square;
xlabel('radius/sqrt(nDim)'); ylabel('density');
title('\bftheory');
legend(legendStrings);


%% simulated densities
n=1e5;
i=1;
stds=[];

subplot(2,2,2);

for nDim=nsDim
    points=randn(n,nDim);
    %points=randsphere(n,nDim,1);
    radii=sqrt(sum(points.^2,2));
    stds(i)=std(radii);
    [counts,xout]=hist(radii,50);
    plot(xout/sqrt(nDim),counts*sqrt(nDim),'Color',cols(i,:),'LineWidth',lw); hold on;
    i=i+1;
end
axis([0 3 0 n/sqrt(2)]); axis square;
xlabel('radius/sqrt(nDim)'); ylabel('count');
title('\bfsimulation');

subplot(2,1,2);
plot(nsDim,stds,'k');
xlabel('number of dimensions');
ylabel('standard deviation of radii');

