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
set(gcf,'Name',sprintf('UniformLikeSum_%uD',d));
xlabel('radius');ylabel('log likelihood')
figure;imagesc(rs,cscales,nansum(est,3))
set(gcf,'Name',sprintf('UniformLike_%uD',d));
xlabel('radius'); ylabel('distance from best center')
title('log likelihood'); colorbar
figure;plot(rs,nansum(est,3)');
set(gcf,'Name',sprintf('UniformLikeMargs_%uD',d));
xlabel('radius');ylabel('log likelihood')

printFig(gcf-2:gcf,[],'pdf')

return


%% DeBUG
nc = 100; nr = 151; n = 200; N = 100;
est = NaN(nc,nr,N,10); rs = linspace(0.5,3,nr);
for d = 2:2:20
   est(:,:,:,d/2) = demo_uniformLikelihood(n,d,N);
   goodr(d/2)=rs(maxix(nansum(nansum(est(:,:,:,d/2),3))));
   [a,b]=ind2sub([size(est,1) size(est,2)],maxix(vectify(nansum(est(:,:,:,d/2),3))));
   bestr(d/2)=rs(b);
end
figure;hold on;plot(2:2:20,bestr);plot(2:2:20,goodr,'g-')
set(gcf,'Name','goodbestrs'); printFig;


