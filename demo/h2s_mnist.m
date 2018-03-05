% h2s_mnist:
% runs hypersphere2sphere on MNIST raw images and trained NN activations
% 
% 2018-03-02 AZ Created

images = loadMNISTImages('data/train-images-idx3-ubyte');
labels = loadMNISTLabels('data/train-labels-idx1-ubyte');
chunk = 1:1000;

label2resp = @(L) cell2mat(arrayfun(@(i) L(chunk)==mod(i,10),1:10,'UniformOutput',false));

categories.vectors = label2resp(labels(chunk));
categories.labels = arrayfun(@num2str,mod(1:10,10),'UniformOutput',false)';
categories.colors = [250   138   117
                     246   136   159
                     215   148   196
                     162   166   215
                     100   180   210
                      61   187   180
                      88   189   139
                     135   184    98
                     182   173    74
                     223   157    79]/255;
titlestr = sprintf('MNIST %u images',numel(chunk));

% MNIST, raw
h=figure(101);clf;
HS2SandMDS(images(:,chunk)',categories,[h 1 4 1],titlestr,2)
h=figure(102);clf;
HS2SandMDS(images(:,chunk)',categories,[h 1 4 1],titlestr,3)

% MNIST, trained model responses, from pyTorch
modelcats = load('data/mnist_fcBNx3.mat'); % Trained 3 layer FC net on MNIST
modelcats.labels = categories.labels;
modelcats.colors = categories.colors;
%modelcats.images = reshape(modelcats.images,[],28*28)';
modelcats.vectors= label2resp(mod(maxix(modelcats.outputs(chunk,:),[],2),10));

if max(modelcats.outputs(:)) <= 0
   modelcats.outputs= modelcats.outputs - min(modelcats.outputs(:)) + 1;
end

h=figure(103);clf;
HS2SandMDS(modelcats.outputs(chunk,:),modelcats,[h 1 4 1],titlestr,2)
h=figure(104);clf;
HS2SandMDS(modelcats.outputs(chunk,:),modelcats,[h 1 4 1],titlestr,3)



