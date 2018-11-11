%demo_histograms

n = 1000000;
npasses = 100;
N = 1;
types = {'Uniform','Gaussian','Hypercubic'};
maxbin = 2;
bins = -0.01:0.01:maxbin+0.01;

fp.DESTINATION = {'elife',[1 2]};
figsetup;fp.txt{4}='Oblique';fp.signif=2;
fh = newfigure('histograms',[1 numel(types)]);

% Dimensions to plot
ds = [1 2 3 4 8 32 128 1024];
% Generate perceptually even light gray -> black ...
cols = flip(makeColorSaturationSeries([0 0 0],numel(ds)));

for i = 1:numel(ds)
   d = ds(i);
   fprintf('Generating histograms for %u dimensions...\n',d)
   for t = 1:numel(types)
      switch t
      case 1,   r = 1;
      case 2,   r = sqrt(2)*exp(gammaln((d+1)/2)-gammaln(d/2));
      case 3,   r = hypercubeRadius(d);
      end
      % Initialize histogram bins, for parallel processing
      h    = zeros(npasses,numel(bins)); 
      parfor j = 1:npasses
         if i==numel(ds), stationarycounter(j,npasses); end
         X      = sampleSpheres(n,d,N,types{t});
         radii  = sqrt(sum(X.^2,2));
         h(j,:) = hist(radii(:)/r,bins);
      end
      % Plot
      h = sum(h);
      h([1 end])=NaN; % don't plot huge tail
      axtivate(t);
      %plot(bins,h/sum(h),'Color',cols(i,:),'LineWidth',2);
      plot(bins(2:end-1),h(2:end-1)/max(h),'Color',cols(i,:),'LineWidth',1);
   end
end

%legend(num2str(vectify([ds;ds])))
%legend(num2str(ds(:)));
axtivate(1);
legend(num2str(ds(:)));
ylabel('Frequency')
axtivate(2)
xlabel('Normalized distance from center')
%arrayfun(@(a) axesSeparate(a,fp),fh.a.h);

batchPlotRefine(fh,fp);
% fp.pap.sz = [3.375 3.375/2];
%set(fh.f,'Renderer','painters','PaperUnits','inches','PaperSize',papsz,...
%    'PaperPosition',[0.01*fp.pap.sz fp.pap.sz],'PaperPositionMode','manual');
delete(separateAxisElements(fh.a.h,'phyplot_y'))
printFig;

