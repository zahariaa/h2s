function h2s_radii(n,d,type)

N = 100; % bootstraps

for d = d
   X  = sampleSpheres(n,d,N,type);
   radii = sqrt(sum(X.^2,2));
   % 0.5*CDF method
   radii = sort(radii,1,'ascend');
   cdf = cumsum(~~radii,1)/n;
   mx = radii(minix(abs(0.5-cdf)))*2^(1/d);
   mx = mean(mx(:));
   md = mean(median(radii,1)   );
%   mx = mean(max(   radii,[],1));
   % expDistPerRad method
   [~,~,mn] = arrayfun(@(i) estimateHypersphere(X(:,:,i)),1:N,'UniformOutput',false);
   mn = mean(cell2mat(mn));

   figure;clf;
   subplot(1,3,2:3);hold on;
   [h,c]=hist(radii(:));bar(c,h);
   plot(squeeze(radii),1.1*max(h)*squeeze(cdf),'k-')
   plot([mn mn],[0 1.1*max(h)],'r-')
   plot([mx mx],[0 1.1*max(h)],'g-')
   plot([md md],[0 1.1*max(h)],'b-')
   xlabel('L2 norm')
   title({type;
          sprintf('%u points in %u dimensions',n,d)});
   subplot(1,3,1);hold on;
   plot(X(:,1,1),X(:,2,1),'ko');
   plot([0 mn]./sqrt(2),[0 mn]./sqrt(2),'r-','LineWidth',8);
   plot([0 mx]./sqrt(2),[0 mx]./sqrt(2),'g-','LineWidth',4);
   plot([0 md]./sqrt(2),[0 md]./sqrt(2),'b-','LineWidth',2);
   axis equal square
   drawCircle;
   title({'r=expDistPreRad method,';
          'g=CDF*2^{1/d}, b=median' })
end
1-[mn md mx]

return

end


