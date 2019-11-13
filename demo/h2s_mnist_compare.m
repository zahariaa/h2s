% h2s_mnist:
% runs hypersphere2sphere on MNIST raw images and trained NN activations
% 
% 2018-03-02 AZ Created

dimLow = 3;
model  = 'lenet';
epsilon = 0.1;
colors = [250   138   117
          246   136   159
          215   148   196
          162   166   215
          100   180   210
           61   187   180
           88   189   139
          135   184    98
          182   173    74
          223   157    79]/255;
% MNIST, trained model responses, from pyTorch
data = load(sprintf('demo/mnist_%s_eps%0.2f.mat',model,epsilon)); % Trained 3 layer FC net on MNIST
label2resp = @(L,chunk) cell2mat_concat(arrayfun(@(i) L(chunk)==i,0:9,'UniformOutput',false))';

chunk = 1:10000;

%% reorganize data structure
nlayers = numel(fieldnames(data))-1;
for i = 1:nlayers
   eval(sprintf('data.x{%u} = data.x%u;data = rmfield(data,"x%u");',i,i,i));
end
%% compute final output
data.x{nlayers+1} = label2resp(maxix(data.x{nlayers}(chunk,:),[],2)'-1,chunk);

%% Construct labels
% right now this uses the true labels
cats.(model) = Categories(label2resp(data.labels,chunk),...
                        arrayfun(@num2str,0:9,'UniformOutput',false)',colors);
% % output labels
% cats.(model) = Categories(data.x{end},...
%                         arrayfun(@num2str,0:9,'UniformOutput',false)',colors);

%% Estimate hyperspheres, run h2s
for i = 1:nlayers+1
   stationarycounter(i,nlayers+1)
   hi.(model)(i) = SetOfHyps('estimate',data.x{i}(chunk,:),cats.(model));
   [hi.(model)(i).sig,hi.sec(i)] = hi.(model)(i).significance(data.x{i}(chunk,:),100);
   lo.(model)(i) = hi.(model)(i).h2s(dimLow);
   if all(lo.(model)(i).radii==0)
      lo.(model)(i).radii = 0.05*ones(1,10);
   end
   [lo.(model)(i).sig,lo.sec(i)] = lo.(model)(i).significance(data.x{i}(chunk,:),100);
end

%% PLOT
subplots = [floor(sqrt(nlayers+1)) ceil(sqrt(nlayers+1))];
fh = newfigure(subplots,[],[],sprintf(       '%s_eps%0.2f_h2s%u',model,epsilon,dimLow));
fh = newfigure(subplots,fh,[],sprintf('values_%s_eps%0.2f_h2s%u',model,epsilon,dimLow));
fh = newfigure(subplots,fh,[],sprintf( 'stats_%s_eps%0.2f_h2s%u',model,epsilon,dimLow));
fh = newfigure(subplots,fh,[],sprintf('stats2_%s_eps%0.2f_h2s%u',model,epsilon,dimLow));
for i = 1:nlayers+1
   lo.(model)(i).show(fh.a(1).h(i));                     title(sprintf('Layer %u',i));
   lo.(model)(i).showValues(fh.a(2).h(i),hi.(model)(i)); title(sprintf('Layer %u',i));
   hi.(model)(i).showSig(fh.a(3).h(i),'legend');         title(sprintf('Layer %u',i));
   hi.(model)(i).showSig(fh.a(4).h(i),hi.sec(i));        title(sprintf('Layer %u',i));
end

if dimLow==2,   printFig(fh,[],'eps')
else            printFig(fh.f(2:4),[],'eps'); printFig(fh.f(1),[],'png',300)
end

