function h2s_radii(n,ds,type)

N = 100; % bootstraps

s = NaN(N,5,numel(ds)); % container for sampled estimators
colors = {[1 0 0],[0 1 0],[0 0 1],[1 1 0],[1 0 1]};

for d = ds
   dix = find(d==ds);
   % Sample & measure radii
   X  = sampleSpheres(n,d,N,type);
   radii = sqrt(sum(X.^2,2));
   maxradii = max(radii,[],1);
   % Estimators
   % expDistPerRad method
   [~,~,tmp] = arrayfun(@(i) estimateHypersphere(X(:,:,i)),1:N,'UniformOutput',false);
   s(:,1,dix) = cell2mat(tmp);
   s(:,3,dix) = median(radii,1);
   s(:,2,dix) = s(:,3,dix).*2^(1/d);   % median*2^(1/d)
   s(:,4,dix) = maxradii + maxradii/(n*(d-1));                      % MVUE for uniform distribution
   v = std(reshape(X,[n*d N])); % squeeze(mean(std(X,[],2))); %v(i) = max(diag(cov(X(:,:,i))));
   s(:,5,dix) = v*sqrt(2)*exp(gammaln((d+1)/2)-gammaln(d/2)); % MVUE for gaussian distribution
   % Results
   figure(98);clf;plotEstimators(X,radii,mat2cell(s(:,:,dix),100,ones(1,5)),colors);
%   skewness(radii(:))
end
1-[mn md mx]

return

end


