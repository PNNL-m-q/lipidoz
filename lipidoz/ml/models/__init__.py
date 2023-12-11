""" 
lipidoz/ml/models/__init__.py

Dylan Ross (dylan.ross@pnnl.gov)

    subpackage for defining individual ML models 
"""


import copy
import time

import torch
from torch import nn, optim
from torch.utils.data import DataLoader
from torchvision import transforms
import numpy as np


class _Model():
    """ 
    Internal base class for dealing with models. Implements the following methods:

    * load -- loads the parameters for a model that has been trained 
    * save -- saves the current model's parameters as a state dict to file
    * train -- train model using specified data and optimization parameters
    * predict -- takes input data and returns predictions

    Subclasses must implement:

    * __init__ -- initializes a new model, untrained. Must set self.model and self.device attributes
    """

    def load(self, state_dict_path):
        """
        Loads parameters for ``self.model`` 

        ``self.model`` has instance of untrained model, is updated with parameters loaded from state dict then gets 
        sent to ``self.device``

        Parameters
        ----------
        state_dict_path : ``str``
            path to state dict with pre-trained parameters for this model
        """
        self.model.load_state_dict(torch.load(state_dict_path))
        self.model = self.model.to(self.device)

    def save(self, state_dict_path):
        """
        Saves parameters for ``self.model`` to file as state dict

        Parameters
        ----------
        state_dict_path : ``str``
            path to state dict to save parameters for this model
        """
        torch.save(self.model.state_dict(), state_dict_path)

    def _train_epoch(self, dataloaders, dataset_sizes, criterion, optimizer, scheduler, 
                     best_loss, best_acc, best_model_wts, debug):
        """
        performs 1 epoch of model training

        Parameters
        ----------
        dataloaders : ``dict(str:torch.utils.data.DataLoader)``
            'train' and 'validate' Dataloaders (``torch.utils.data.DataLoader``)
        dataset_sizes : ``dict(str:int)``
            'train' and 'validate' dataset sizes
        criterion : ``torch.nn.?``
            loss function for training
        optimizer : ``torch.optim.Optimizer``
            model optimizer
        scheduler : ``torch.optim.lr_schedulerStepLR``
            learning rate scheduler
        best_loss : ``float``
            current minimum loss value
        best_acc : ``float``
            current maximum accuracy value
        best_model_wts : ``?``
            current best performing model weights
        debug : ``bool``
            print debugging info

        Returns
        -------
        best_loss : ``float``
            new minimum loss value
        best_acc : ``float``
            new maximum accuracy value
        best_model_wts : ``?``
            new best performing model weights
        """
        if debug:
            # keep track of elapsed time in epoch
            t0 = time.time()
            print('learning rate: {:.5f}'.format(scheduler.get_last_lr()[0]))
        # Each epoch has a training and validation phase
        for phase in ['train', 'val']:
            if phase == 'train':
                self.model.train()  # Set model to training mode
            else:
                self.model.eval()   # Set model to evaluate mode
            # keep track of loss/corrects
            running_loss = 0.0
            running_corrects = 0
            # Iterate over data.
            for inputs, labels in dataloaders[phase]:
                inputs = inputs.to(self.device)
                labels = labels.to(self.device)
                # zero the parameter gradients
                optimizer.zero_grad()
                # forward
                # track history if only in train
                with torch.set_grad_enabled(phase == 'train'):
                    outputs = self.model(inputs.float())
                    _, preds = torch.max(outputs, 1)
                    loss = criterion(outputs, labels)
                    # backward + optimize only if in training phase
                    if phase == 'train':
                        loss.backward()
                        optimizer.step()
                # statistics
                running_loss += loss.item() * inputs.size(0)
                running_corrects += torch.sum(preds == labels.data)
            if phase == 'train':
                scheduler.step()
            epoch_loss = running_loss / dataset_sizes[phase]
            epoch_acc = running_corrects.double() / dataset_sizes[phase]
            if debug:
                if phase == 'train':
                    msg = '  phase: {}\n      epoch loss: {:.2e}\n      epoch acc:  {:.4f}'
                    print(msg.format(phase, epoch_loss, epoch_acc))
                else:
                    msg = '  phase: {}\n      epoch loss: {:.2e} (best: {:.2e})\n      epoch acc:  {:.4f}   (best: {:.4f})'
                    print(msg.format(phase, epoch_loss, best_loss, epoch_acc, best_acc))
            if phase == 'val' and epoch_loss < best_loss:
                best_loss = epoch_loss
                best_acc = epoch_acc
                best_model_wts = copy.deepcopy(self.model.state_dict())
        if debug:
            # report elapsed time in epoch
            elapsed = time.time() - t0
            print('epoch time: {:.0f}m {:.0f}s'.format(elapsed // 60, elapsed % 60))
        return best_loss, best_acc, best_model_wts

    def train(self, dataloaders, dataset_sizes, 
              criterion=None, optimizer=None, scheduler=None, epochs=32, debug=False, xent_f_t_weights=[0.9, 0.1]):
        """
        Trains a model with specified dataset and training parameters

        Parameters
        ----------
        dataloaders : ``dict(str:torch.utils.data.DataLoader)``
            'train' and 'validate' Dataloaders (``torch.utils.data.DataLoader``)
        dataset_sizes : ``dict(str:int)``
            'train' and 'validate' dataset sizes
        criterion : ``torch.nn.?``, optional
            loss function for training, if not provided defaults to ``torch.nn.CrossEntropyLoss`` with weight of 
            True examples set to 10% and False to 90% to reflect the imbalance in the training data
        optimizer : ``torch.optim.Optimizer``, optional
            model optimizer, if not provided defaults to ``torch.optim.Adam`` with learning rate of 0.001
        scheduler : ``torch.optim.lr_schedulerStepLR``, optional
            learning rate scheduler, if not provided defaults to decaying learning rate by 0.1 every 8 epochs
        epochs : ``int``, default=32
            number of epochs to train over
        debug : ``bool``, default=False
            print debugging info
        xent_f_t_weights : ``list(float)``, default=[0.9, 0.1]
            if using the default cross-entropy loss, set weights for [F, T] classes. By defalt this
            ratio is 0.9 F to 0.1 T to reflect the approximate imbalance in training examples, but 
            the ratio can be tuned to achieve desired prediction characteristics
        """
        if debug:
            # keep track of elapsed time
            t0 = time.time()
            print('training model...')
            msg = 'training examples: {train} validation examples: {val}'
            print(msg.format(**dataset_sizes))
        # deal with all of the training utilities
        if criterion is None:
            # weight False/True 90%/10% to approximate imbalance of labels in training data
            # this is now configurable
            criterion = torch.nn.CrossEntropyLoss(weight=torch.tensor(xent_f_t_weights))
        if optimizer is None:
            # Adam optimizer, 0.001 learning rate
            optimizer = torch.optim.Adam(self.model.parameters(), lr=0.001)
        if scheduler is None:
            # decay learning rate by 0.5 every 4 epochs
            scheduler = optim.lr_scheduler.StepLR(optimizer, step_size=4, gamma=0.5)
        # proceed with training
        best_model_wts = copy.deepcopy(self.model.state_dict())
        best_acc = 0.0
        best_loss = 1e9
        for epoch in range(epochs):
            if debug:
                print('---------- epoch ({:2d}/{:2d}) ----------'.format(epoch + 1, epochs))
            best_loss, best_acc, best_model_wts = self._train_epoch(dataloaders, dataset_sizes, 
                                                                    criterion, optimizer, scheduler, 
                                                                    best_loss, best_acc, best_model_wts, debug)
        if debug:
            print('------------------------------------')
            elapsed = time.time() - t0
            print('total training time: {:.0f}m {:.0f}s'.format(elapsed // 60, elapsed % 60))
            print('best val acc:  {:.4f}'.format(best_acc))
            print('best val loss: {:.2e}'.format(best_loss))
        # load best model weights, set self.model
        self.model.load_state_dict(best_model_wts)

    def _get_dataloader_from_input(self, X):
        """ takes raw data array and makes a dataloader for it, fixes shape if it is in pre-ml format """
        X_ = []
        for mat in X:
            if mat.shape[-1] != 3:
            # account for data that is shaped like pre-ml data 
                mat = mat.transpose(1, 2, 0)
            mat = transforms.Compose([transforms.ToTensor()])(mat)
            X_.append(mat)
        return DataLoader(X_, batch_size=256, shuffle=False)

    def predict(self, X):
        """
        Predict class labels for input examples

        Parameters
        ----------
        X : ``numpy.ndarray``
            array of input data for ML with shape (N, 24, 400, 3), where N is the number of examples in the set. 
            Shape (N, 3, 24, 400) (as in pre-ml data) is also ok, it automatically gets transposed to the proper shape

        Returns
        -------
        y : ``numpy.ndarray``
            array of predictions, 0 for False 1 for True, with shape (N,) where N is the number of examples in the set
        """
        # make a dataloader from the input data
        dataloader = self._get_dataloader_from_input(X)
        was_training = self.model.training
        y = np.array([])
        self.model.eval()
        with torch.no_grad():
            for inputs in dataloader:
                inputs = inputs.to(self.device)
                outputs = self.model(inputs.float())
                _, preds = torch.max(outputs, 1)
                y = np.concatenate([y, preds])
        self.model.train(mode=was_training)
        return y

    def predict_proba(self, X):
        """
        Predicts class probabilities for input examples
        
        Parameters
        ----------
        X : ``numpy.ndarray``
            array of input data for ML with shape (N, 24, 400, 3), where N is the number of examples in the set. 
            Shape (N, 3, 24, 400) (as in pre-ml data) is also ok, it automatically gets transposed to the proper shape

        Returns
        -------
        y : ``numpy.ndarray``
            array of label (T/F) probabilities, with shape (N, 2) where N is the number of examples in the set
        """
        # make a dataloader from the input data
        dataloader = self._get_dataloader_from_input(X)
        was_training = self.model.training
        y = []
        self.model.eval()
        with torch.no_grad():
            for inputs in dataloader:
                inputs = inputs.to(self.device)
                outputs = self.model(inputs.float())
                probs = nn.functional.softmax(outputs, dim=1)
                y += probs.tolist()
        self.model.train(mode=was_training)
        return np.array(y)
