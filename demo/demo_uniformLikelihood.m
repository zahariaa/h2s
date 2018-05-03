function est = demo_uniformLikelihood(n,d,N)

if ~exist('n','var') || isempty(n),   n = 200;   end
if ~exist('d','var') || isempty(d),   d =  20;   end
if ~exist('N','var') || isempty(N),   N = 100;   end

X=sampleSpheres(n,d,1,'Uniform');

nr = 151; nc = 100;
rs = linspace(0.5,3,nr);
cscales = linspace(0,2,nc);
est = NaN(nc,nr,N);

%% Bootstrap
for b = 1:N
   stationarycounter(b,N);
   cnoise = randn(N,d);
   cnoise = cnoise./repmat(sqrt(sum(cnoise.^2,2)),[1 d]);
   cnoise = cnoise.*repmat(cscales',[1 d]);
   cs = repmat(mean(X),[nc 1]) + cnoise;
   for i = 1:nc
      for j = 1:nr
est(i,j,b) = log(uniformLikelihood([rs(j) cs(i,:)],X));
      end
   end
end
est(isinf(est)) = NaN;
est(est == 0  ) = NaN;

%% Plot
figure;plot(rs,nansum(nansum(est,3)));
set(gcf,'Name',sprintf('%uD-UniformLikeSum',d));
xlabel('radius');ylabel('log likelihood')
figure;imagesc(rs,cscales,nansum(est,3))
set(gcf,'Name',sprintf('%uD-UniformLike',d));
xlabel('radius'); ylabel('distance from best center')
title('log likelihood'); colorbar
figure;plot(rs,nansum(est,3)');
set(gcf,'Name',sprintf('%uD-UniformLikeMargs',d));
xlabel('radius');ylabel('log likelihood')

printFig(gcf-2:gcf,[],'pdf')

return


%% DeBUG
nc = 100; nr = 151; n = 200; N = 100;
est = NaN(nc,nr,N,10);
for d = 2:2:20
   est(:,:,:,d/2) = demo_uniformLikelihood(n,d,N);
   [a,b]=ind2sub([size(est,1) size(est,2)],maxix(nansum(nansum(est(:,:,:,d/2),3))));
   goodr(d/2)=rs(b);
end

