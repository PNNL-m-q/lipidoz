""" 
lipidoz/ml/models/resnet18.py

Dylan Ross (dylan.ross@pnnl.gov)

    deep-learning model based on pre-trained RESNET18
"""


import torch
from torch import nn
from torchvision import models

from lipidoz.ml.models import _Model


class ResNet18(_Model):
    """
    Model based on pre-trained ResNet18
    """

    def __init__(self):
        """
        Inits a new instance of RESNET18 model
        """
        # use CUDA device if available otherwise use the CPU
        self.device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
        # set up the model
        self.model = models.resnet18(pretrained=True)
        num_ftrs = self.model.fc.in_features
        self.model.fc = nn.Linear(num_ftrs, 2)
        self.model = self.model.to(self.device)

