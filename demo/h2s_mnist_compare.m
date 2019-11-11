% h2s_mnist:
% runs hypersphere2sphere on MNIST raw images and trained NN activations
% 
% 2018-03-02 AZ Created

dimLow = 3;
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
lenet = load('demo/mnist_LeNet5.mat'); % Trained 3 layer FC net on MNIST
label2resp = @(L,chunk) cell2mat_concat(arrayfun(@(i) L(chunk)==i,0:9,'UniformOutput',false))';

chunk = 1:10000;

%% reorganize data structure
for i = 1:5
   eval(sprintf('lenet.x{%u} = lenet.x%u;lenet = rmfield(lenet,"x%u");',i,i,i));
end
%% compute final output
lenet.x{6} = label2resp(maxix(lenet.x{5}(chunk,:),[],2)'-1,chunk);

%% Construct labels
% right now this uses the true labels
cats.lenet = Categories(label2resp(lenet.labels,chunk),...
                        arrayfun(@num2str,0:9,'UniformOutput',false)',colors);
% % output labels
% cats.lenet = Categories(lenet.x{6},...
%                         arrayfun(@num2str,0:9,'UniformOutput',false)',colors);

%% Estimate hyperspheres, run h2s
for i = 1:6
   hi.lenet(i) = SetOfHyps('estimate',lenet.x{i}(chunk,:),cats.lenet);
   lo.lenet(i) = hi.lenet(i).h2s(dimLow);
   if all(lo.lenet(i).radii==0)
      lo.lenet(i).radii = 0.05*ones(1,10);
   end
   [lo.lenet(i).sig,lo.sec(i)] = lo.lenet(i).significance(lenet.x{i}(chunk,:),100);
end

%% PLOT
fh = newfigure(sprintf('lenet_h2s%u',dimLow),[2 3]);
for i = 1:6; axtivate(fh.a.h(i)); lo.lenet(i).show; end

if dimLow==2,   printFig(fh.f,[],'eps')
else            printFig(fh.f,[],'png',450)
end

