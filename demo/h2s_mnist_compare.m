function h2s_mnist_compare(epsilon)
% h2s_mnist:
% runs hypersphere2sphere on MNIST raw images and trained NN activations
% 
% 2018-03-02 AZ Created

dimLow = 2;
model  = 'lenet';
%epsilon = 0;
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
h2sfile = sprintf('demo/mnist_%s_eps%0.2f_h2s.mat',model,epsilon);
if exist(h2sfile,'file'), load(h2sfile); fprintf('Loaded h2s from %s\n',h2sfile);
else
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
   [hi.(model)(i).sig,hi.sec(i)] = hi.(model)(i).significance(data.x{i}(chunk,:),1000);
   lo.(model)(i) = hi.(model)(i).h2s(dimLow);
   if all(lo.(model)(i).radii==0)
      lo.(model)(i).radii = 0.05*ones(1,10);
   end
   [lo.(model)(i).sig,lo.sec(i)] = lo.(model)(i).significance(data.x{i}(chunk,:),1000);
end
save(h2sfile,'hi','lo','cats','nlayers')
end

%% PLOT
fh = newfigure([3 nlayers+1],[],[],sprintf( 'combo_%s_eps%0.2f_h2s%u',model,epsilon,dimLow));
fh = newfigure([1 nlayers+1],fh,[],sprintf('stats2_%s_eps%0.2f_h2s%u',model,epsilon,dimLow));
for i = 1:nlayers+1
   if i<nlayers+1,   titlestr = sprintf('Layer %u',i);
   else              titlestr = 'Final readout';
   end
   lo.(model)(i).show(fh.a(1).h(i));                               title(titlestr);
   lo.(model)(i).showValues(fh.a(1).h(i+nlayers+1),hi.(model)(i));
   hi.(model)(i).showSig(fh.a(1).h(i+2*(nlayers+1)));
   hi.(model)(i).showSig(fh.a(2).h(i),hi.sec(i));                  title(titlestr);
end

if dimLow==2,   printFig(fh.f,[],'eps')
else            printFig(fh.f(2:end),[],'eps'); printFig(fh.f(1),[],'png',300)
end

