function L = uniformLikelihood(rad_center,X)
% rad_center: first element is radius, rest is center

%% Preliminaries
rad    = rad_center(1);
center = rad_center(2:end);
[n,d]  = size(X);

%% Calculations
volsphere = ((rad^d)*(pi^(d/2)))/gamma((d+1)/2);
radii = sqrt(sum((X - repmat(center(:)',[n 1])).^2,2));

L = sum(radii<=rad)/volsphere;

return


%% DEBUG
n=200;d=8;N=1;X=sampleSpheres(n,d,N,'Uniform');
rs = 0.5:0.01:2;
for b = 1:100 %bootstrap 
   stationarycounter(b,100);
   cnoise = randn(100,d);
   cnoise = cnoise./repmat(sqrt(sum(cnoise.^2,2)),[1 d]);
   cnoise = cnoise.*repmat(linspace(0,2,100)',[1 d]);
   cs = repmat(mean(X),[100 1]) + cnoise;
   for i = 1:100
      for j = 1:numel(rs)
est(i,j,b) = log(uniformLikelihood([rs(j) cs(i,:)],X));
      end
   end
end
est(isinf(est))=NaN;
figure;plot(rs,nansum(nansum(est,3)));
set(gcf,'Name',sprintf('%uD-UniformLikeSum',d));
xlabel('radius');ylabel('log likelihood')
figure;imagesc(rs,linspace(0,2,100),nansum(est,3))
set(gcf,'Name',sprintf('%uD-UniformLike',d));
xlabel('radius'); ylabel('distance from best center')
title('log likelihood'); colorbar


