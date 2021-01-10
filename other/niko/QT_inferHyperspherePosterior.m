% QT_inferHyperspherePosterior

nSamples = 10000;

%% 1-d scenario
%points = 0;
%points = [-1 1]';
points = linspace(-2,2,5)';

%points = [-9 -8 -7 9]';

[posteriorSamplesLoc,posteriorSamplesRad,logLikelihoods] = inferHyperspherePosterior(points, nSamples);


%% Better visualization
c=posteriorSamplesLoc; r=posteriorSamplesRad;
L=logLikelihoods; L=(L-min(L))/(max(L)-min(L));

% Threshold likelihood color
thresh = 0.5; L = (L-thresh)/(1-thresh); L(L<0) = 0;

[L,ix]=sort(L);
r = r(ix); c = c(ix);
ix = round(linspace(nSamples/9,nSamples,9));

fh = newfigure([2 2],'jointPosterior');
axtivate(1);
plot(points,0,'wo','MarkerFaceColor','k');
for i = 1:2:numel(ix)-1
   plot(c(ix(i))+[-1 1]*r(ix(i)),max(ix)*1.1-[ix(i) ix(i)],...
        '-','Color',L(ix(i))*[1 0 0],'LineWidth',8)
end
plot(median(c)+[-1 1]*median(r),max(ix)*[1.1 1.1],...
     '-b','LineWidth',8)
plot(c(end)   +[-1 1]*r(end)   ,max(ix)*[0.1 0.1],...
     '-g','LineWidth',8)
box off; fh.a.h(1).YAxis.Visible='off';
pos = get(fh.a.h(1),'Position');
set(fh.a.h(1),'Position',pos.*[1 1 0.5 0.5]);

axtivate(4)
for i=1:numel(L)
   plot(c(i),r(i),'o',...
   'MarkerFaceColor',L(i)*[1 0 0],'MarkerEdgeColor','none');
end
plot(median(c),median(r),'wo',...
     'MarkerFaceColor','b','MarkerSize',8)
plot(c(end),r(end),'wo',...
     'MarkerFaceColor','g','MarkerSize',8)
for i = 1:2:numel(ix)-1
   plot(c(ix(i)),r(ix(i)),'o','MarkerFaceColor','w',...
       'MarkerEdgeColor',L(ix(i))*[1 0 0],'MarkerSize',8)
end
box off
xlabel('Center estimate')
ylabel('Radius estimate')
set(fh.a.h(4),'YAxisLocation','right')

%% TODO: colorbar with thresh
%colorbar

axtivate(3);
redges = linspace(floor(min(r)),ceil(max(r)),40);
hr = histogram(r,redges,'FaceColor','k','EdgeColor','w','FaceAlpha',1);
hr.Orientation = 'horizontal';
fh.a.h(3).YAxis.Visible = 'off';
fh.a.h(3).XDir          = 'reverse';
% Convert counts to frequency
fh.a.h(3).XTickLabel = num2cell(num2str(cellfun(@str2double,...
                       fh.a.h(3).XTickLabel)/nSamples),2);
xlabel('Frequency')
pos = get(fh.a.h(3),'Position');
set(fh.a.h(3),'Position',pos.*[1 1 0.5 1]);
ylim(ylim(fh.a.h(4)));

axtivate(2);
cedges = linspace(floor(min(c)),ceil(max(c)),40);
hc = histogram(c,cedges,'FaceColor','k','EdgeColor','w','FaceAlpha',1);
fh.a.h(2).XAxis.Visible = 'off';
fh.a.h(2).YTickLabel = num2cell(num2str(cellfun(@str2double,...
                       fh.a.h(2).YTickLabel)/nSamples),2);
ylabel('Frequency')
pos = get(fh.a.h(2),'Position');
set(fh.a.h(2),'Position',pos.*[1 1 1 0.5]);
set(fh.a.h(2),'YAxisLocation','right')
xlim(xlim(fh.a.h(4)));

printFig([],[],'eps')

