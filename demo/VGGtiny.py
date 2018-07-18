"""
Reduced VGG-16 for Tiny ImageNet Model
modified VGG-16, code mostly copied from torchvision.model
"""
import math
import torch
import torch.nn as nn
import torch.optim as optim
# For dynamic loss plots
import matplotlib
import matplotlib.pyplot as plt
from IPython.display import clear_output

class VGGtiny(nn.Module):

    def __init__(self, features, num_classes=200, init_weights=True):
        super(VGGtiny, self).__init__()
        self.features = features
        self.classifier = nn.Sequential(
            nn.Linear(512 * 7 * 7, 4096),
            nn.ReLU(True),
#             nn.Dropout(),
            nn.Linear(4096, 2048),         # reduced from 4096,4096
            nn.ReLU(True),
#             nn.Dropout(),
            nn.Linear(2048, num_classes),  # reduced from 4096, num_classes
        )
        if init_weights:
            self._initialize_weights()

    def forward(self, x):
        x = self.features(x)
        x = x.view(x.size(0), -1)
        x = self.classifier(x)
        return x

    def _initialize_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Conv2d):
                nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0)
            elif isinstance(m, nn.BatchNorm2d):
                nn.init.constant_(m.weight, 1)
                nn.init.constant_(m.bias, 0)
            elif isinstance(m, nn.Linear):
                nn.init.normal_(m.weight, 0, 0.01)
                nn.init.constant_(m.bias, 0)

    def train(self, trainloader, plotmode=False, nupdate=5):
        # zero gradients
        optimizer = optim.SGD(self.parameters(), lr=0.01, momentum=0.9)
        optimizer.zero_grad()   # zero the gradient buffers

        # set up loss function
        criterion = nn.CrossEntropyLoss()

        for epoch in range(2):  # loop over the dataset multiple times

            running_loss = 0.0
            losses = np.array([])
            epoch_correct = 0.0
            for i, data in enumerate(trainloader, 0):
                # get the inputs
                inputs, labels = data

                # zero the parameter gradients
                optimizer.zero_grad()

                # forward + backward + optimize
                outputs = self(inputs)
                loss = criterion(outputs, labels)
                losses = np.append(losses, loss.item())
                loss.backward()
                optimizer.step()
                
                # calculate accuracy
                epoch_correct += (outputs.max(1)[1] == labels).sum().data[0]

                # print statistics
                running_loss += loss.item()
                if plotmode and i % nupdate == nupdate-1:    # print every 100 mini-batches
        #             print('[%d, %5d] loss: %.3f' %
        #                   (epoch + 1, i + 1, running_loss / 5))
        #             running_loss = 0.0
                    # PLOT!
                    clear_output(wait=True)
                    plt.figure(figsize=(5,5))
                    plt.subplot(1, 2, 1)
                    plt.plot(losses / trainloader.batch_size)
                    plt.xlabel('minibatches')
                    plt.title('loss')
                    
                    plt.subplot(1, 2, 2)
                    plt.plot(epoch_correct / trainloader.batch_size)
                    plt.xlabel('minibatches')
                    plt.title('accuracy')
                    
                    plt.tight_layout()
                    plt.show()

        print('Finished Training')

        

def make_layers(cfg, batch_norm=False):
    layers = []
    in_channels = 3
    for v in cfg:
        if v == 'M':
            layers += [nn.MaxPool2d(kernel_size=2, stride=2)]
        else:
            conv2d = nn.Conv2d(in_channels, v, kernel_size=3, padding=1)
            if batch_norm:
                layers += [conv2d, nn.BatchNorm2d(v), nn.ReLU(inplace=True)]
            else:
                layers += [conv2d, nn.ReLU(inplace=True)]
            in_channels = v
    return nn.Sequential(*layers)

cfg = {
    'A': [64, 'M', 128, 'M', 256, 256, 'M', 512, 512, 'M', 512, 512, 'M'],
    'B': [64, 64, 'M', 128, 128, 'M', 256, 256, 'M', 512, 512, 'M', 512, 512, 'M'],
    'D': [64, 64, 'M', 128, 128, 'M', 256, 256, 256, 'M', 512, 512, 512], # nixed last set: [..., 'M', 512, 512, 512, 'M']
    'E': [64, 64, 'M', 128, 128, 'M', 256, 256, 256, 256, 'M', 512, 512, 512, 512, 'M', 512, 512, 512, 512, 'M'],
}

def vgg16():
    """VGG 16-layer model (configuration "D")
    """
    model = VGGtiny(make_layers(cfg['D']))
    return model

