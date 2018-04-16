function h2s_radii(n,ds,type)

N = 100; % bootstraps

s = NaN(N,5,numel(ds)); % container for sampled estimators
colors = {[1 0 0],[0 1 0],[0 0 1],[1 1 0],[1 0 1]};
target = NaN(1,numel(ds));

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

   m(dix,:) = 1-mean(s(:,:,dix),1);
   if     strcmpi(type,'gaussian'),   target(dix) = median(radii(:));
   elseif strcmpi(type,'uniform' ),   target(dix) = 1;
   end
   if dix > 1
      figure(99);cla;hold on;
      for e = 1:5
         err = (repmat(target,[N 1])-squeeze(s(:,e,:)))./repmat(target,[N 1]);
         plotErrorPatch(log2(ds),log10(err.^2),colors{e});
      end
   end
   xlabel('dimensions')
   ylabel('log_{10} error')
   title([type ' distribution, radius estimation'])
   legend('','expDistPerRad','','median2^{1/d}',...
          '','median','','MVUE-Unif','','MVUE-Gauss',...
          'Location','EastOutside')
end

return

%% DEMOS/DEBUG
h2s_radii(200,log2space(1,12,12),'Uniform');xlim([1 12]);
set(gca,'XTick',1:3:12,'XTickLabel',num2str(2.^(1:3:12)'))
set(99,'Name','RadiusUniform','renderer','painters');
axesSeparate;printFig(99,'~/Desktop/','eps');

h2s_radii(200,log2space(1,12,12),'Gaussian');xlim([1 12]);
set(gca,'XTick',1:3:12,'XTickLabel',num2str(2.^(1:3:12)'))
set(99,'Name','RadiusGaussian','renderer','painters');
axesSeparate;printFig(99,'~/Desktop/','eps');

end


