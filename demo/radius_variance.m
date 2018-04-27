function radius_variance
% radius_variance()
% 
% Demo comparing variance of measured radii as a function
% of dimensionality, for uniform n-ball and gaussian distributions
% 
% 2018-04-19 AZ Created

n = 20;
N = 100;
types = {'Uniform','Gaussian'};

rv = NaN(2,9,N);

for i = 1:9
   d = 2^i;
   targetg = sqrt(2)*exp(gammaln((d+1)/2)-gammaln(d/2));
   for t = 1:2
      switch t
         case 1,   target = 1;
         case 2,   target = targetg;
      end
      X = sampleSpheres(n,d,N,types{t});
      radii = squeeze(sqrt(sum(X.^2,2)));
      rv(t,i,:) = var(radii/target);
   end
end


figure;set(gcf,'Name','radiusvariance');hold on;
plotErrorPatch(1:9,log10(squeeze(rv(1,:,:)))',[0 0 0])
plotErrorPatch(1:9,log10(squeeze(rv(2,:,:)))',[1 0 0])
legend('','Uniform','','Gaussian','Location','SouthWest')
set(legend,'Box','off');
delete(findobj(get(legend,'Children'),'type','patch'));
xlabel('dimensions')
ylabel('log_{10} variance of radii')
set(gca,'XLim',[1 9],'XTick',1:3:9,'XTickLabel',num2str(2.^(1:3:9)'))
axesSeparate;   set(gcf,'Renderer','painters');
printFig([],[],'eps')

