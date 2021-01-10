function showSpheresModelWithDiagnostics(model, radii, distances, margins, figI, titleStr)

%% control variables
patchDetail = 100;
opacity = 1;


%% show spheres model
showSpheresModel(model, [figI 2 2 1], titleStr)


%% summary stats for data and model as matrices
subplot(2,2,2);
showSummaryStats(radii, distances, margins, model)

% summary stats for data and model as bar graphs
dataCol =  [0 0 0];
modelCol = [0.5 0.5 0.5];

subplot(2,3,4);
radiusBars = [radii' model.radii'];
h=bar(radiusBars,'grouped','EdgeColor','none');
colormap([dataCol; modelCol]);
title({'\bfradii',['\rmerror: ',num2str(model.radiiError,2)]});
legend({'data','model'},'Location','SouthEast');

subplot(4,3,[8 9]);
distanceBars = [distances' model.distances'];
if numel(distances)==1, distanceBars = [distanceBars; [0 0]]; end
h=bar(distanceBars,'group','EdgeColor','none');
colormap([dataCol; modelCol]);
title({'\bfdistances',['\rmerror: ',num2str(model.distancesError,2)]});
%legend({'data','model'},'Location','SouthEast');
nanaxis([0.5 numel(distances)+0.5 nan nan]);

subplot(4,3,[11 12]);
marginBars = [margins' model.margins'];
if numel(margins)==1, marginBars = [marginBars; [0 0]]; end
h=bar(marginBars,'group','EdgeColor','none');
colormap([dataCol; modelCol]);
title({'\bfmargins',['\rmerror: ',num2str(model.marginsError,2)]});
%legend({'data','model'},'Location','SouthEast');
nanaxis([0.5 numel(margins)+0.5 nan nan]);


%% show summary statistics
function showSummaryStats(radii, distances, margins, model)

% marginCols = colourScale([ 22 31 40; 255 100 100]/255,100);
% radiiepCols = colourScale([240 240 240; 22 31 40]/255,100);
marginCols = colourScale([ 1 1 1; 1 0 0],100);
radiiepCols = colourScale([0 0 0; 1 1 1],100);

nCats = model.nCats;
marginWidth = nCats^(1/2);
axis([0 nCats 0 (2+marginWidth)*nCats]); axis equal;
set(gca,'YDir','reverse');

margins = squareform(margins);
distances = squareform(distances);

maxradiiep = max([radii(:);model.radii(:);distances(:);model.distances(:)]);
maxmargin = max([margins(:);model.margins(:)]);

shift = 0;
w = 0.9;
lw = 2;
for displayI = 1:2
    for cat1i = 1:nCats
        for cat2i = cat1i:nCats
            if cat1i==cat2i
                % colour patch for radius
                col = radiiepCols(1+round(radii(cat1i)/maxradiiep*99),:);
                rectangle('Position',[cat1i-1,shift+cat2i-1,w,w],'Curvature',0.2,'FaceColor',col,'EdgeColor',model.catCols(cat1i,:),'LineWidth',lw);
            else
                % colour patch for distance
                if maxradiiep == 0
                    col = marginCols(round(size(radiiepCols,1)/2),:);
                else
                    col = radiiepCols(1+round(distances(cat1i,cat2i)/maxradiiep*99),:);
                end
                rectangle('Position',[cat1i-1,shift+cat2i-1,w,w],'Curvature',0.2,'FaceColor',col,'EdgeColor','none');

                % colour patch for margin
                if maxmargin == 0
                    col = marginCols(round(size(marginCols,1)/2),:);
                else
                    col = marginCols(1+round(margins(cat1i,cat2i)/maxmargin*99),:);
                end
                rectangle('Position',[cat2i-1,shift+cat1i-1,w,w],'Curvature',0.2,'FaceColor',col,'EdgeColor','none');
            end
        end % cat2i
    end % cat1i
    
    shift = nCats*(1+marginWidth);
    radii = model.radii;
    distances = squareform(model.distances);
    margins = squareform(model.margins);
    
end % displayI

% labels
txtMarg = 0.3;
text(0, 0-2*txtMarg, '\bfdata');
text(0, shift-2*txtMarg, '\bfmodel');

for catI = 1:nCats
    text(nCats+0.3, catI-1+0.5, model.catLabels{catI});
    text(nCats+txtMarg, shift+catI-1+0.5, model.catLabels{catI});
    text(catI-1+0.5, shift+nCats+txtMarg, model.catLabels{catI},'Rotation',90,'HorizontalAlignment','right');      
end


% colour scales
hold on;
h = 1.0*txtMarg; % height of horizontal colourbars
s = 3*txtMarg; % vertical offset between label and colourbar
s2 = 5*txtMarg; % vertical offset between colourbars

y = shift+nCats+18*txtMarg; % vertical starting point
text(0, y, '\bfradii and distances \rm[min-max]');
image([0 nCats],[y+s y+s+h],reshape(radiiepCols,[1 size(radiiepCols)]));

y = y+s+h+s2;
text(0, y, '\bfmargins \rm[min-max]');
image([0 nCats],[y+s y+s+h],reshape(marginCols,[1 size(marginCols)]));

y = y+s+h+s2;
axis([0 nCats*1.5 -nCats-4 y]);
axis equal off;




