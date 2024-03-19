import numpy as np
import torch
from torch.autograd import Variable
import matplotlib.pyplot as plt
import matplotlib
from sklearn import preprocessing
from sklearn.model_selection import train_test_split
from torch.utils.data import Dataset
from scipy.io import mmread
import math
import sys
from typing import List, Tuple

matplotlib.rcParams.update({'font.size': 17})
plt.style.use("seaborn-whitegrid")
torch.manual_seed(0)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
np.random.seed(0)

device = torch.device('cpu')
batch_size = 5 #8  #10 #64
hidden_size = 3 ###512 ### 256 #500
num_layers =  1  #2
train_size = 0.7 ##0.7

class NutDataset(Dataset):
    def __init__(self, X, y):
        self.X = torch.from_numpy(X).type(torch.float32)
        self.y = torch.from_numpy(y).type(torch.float32)


    def __len__(self):
        return len(self.X)


    def __getitem__(self, idx):
        return [self.X[idx], self.y[idx]]


class lstm_encoder(torch.nn.Module):
    ''' Encodes time-series sequence '''
    def __init__(self, input_size, hidden_size):
       '''
       : param input_size:     the number of features in the input X
       : param hidden_size:    the number of features in the hidden state h
       : param num_layers:     number of recurrent layers (i.e., 2 means there are  2 stacked LSTMs)
       '''
       super(lstm_encoder, self).__init__()
       self.input_size = input_size
       self.hidden_size = hidden_size

        # define LSTM layer
       self.lstm = torch.nn.LSTM(input_size, hidden_size, num_layers=1)

    def forward(self, x_input):
      '''
      : param x_input:               input of shape (seq_len, # in batch, input_size)
      : return lstm_out, hidden:     lstm_out gives all the hidden states in the sequence;
      :                              hidden gives the hidden state and cell state for the last element in the sequence
      '''
      lstm_out, hidden = self.lstm(x_input)
      return lstm_out, hidden           


class lstm_decoder(torch.nn.Module):
    ''' Decodes hidden state output by encoder '''
    def __init__(self, input_size, hidden_size):
        '''
        : param input_size:     the number of features in the input X
        : param hidden_size:    the number of features in the hidden state h
        : param num_layers:     number of recurrent layers (i.e., 2 means there are 2 stacked LSTMs)
        '''
        super(lstm_decoder, self).__init__()
        self.input_size = input_size
        self.hidden_size = hidden_size
        self.lstm = torch.nn.LSTM(input_size, hidden_size, num_layers=1)
        self.linear = torch.nn.Linear(hidden_size, input_size)

    def forward(self, x_input: torch.Tensor, encoder_hidden_states: Tuple[torch.Tensor, torch.Tensor]) -> Tuple[torch.Tensor, Tuple[torch.Tensor, torch.Tensor]]:
        #def forward(self, x_input: torch.Tensor, encoder_hidden_states: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        '''
        : param x_input:                    should be 2D (batch_size, input_size)
        : param encoder_hidden_states:      hidden states
        : return output, hidden:            output gives all the hidden states in the sequence;
        :                                   hidden gives the hidden state and cell state for the last element in the sequence
        '''
        lstm_out, hidden = self.lstm(x_input, encoder_hidden_states)
        output = self.linear(lstm_out[-1, :, :])
        return output, hidden


class lstm_seq2seq(torch.nn.Module):
    ''' train LSTM encoder-decoder and make predictions '''
    def __init__(self, input_size, hidden_size):
     
        super(lstm_seq2seq, self).__init__()
        self.input_size = input_size
        self.hidden_size = hidden_size
        self.encoder = lstm_encoder(input_size, hidden_size)
        self.decoder = lstm_decoder(input_size, hidden_size)

    def forward(self, x):
        # encoder outputs
        encoder_output, encoder_hidden_states = self.encoder(x)
        # x[-1, :, :]   # shape: (batch_size, input_size) because (#seq_len, #batch_size, #input_size) if (#batch_size, #seq_len, #input_size)
        decoder_input = x
        decoder_hidden = encoder_hidden_states
        decoder_output, decoder_hidden = self.decoder(decoder_input, decoder_hidden)
        decoder_input = decoder_output
        return decoder_output
class Net:
    def __init__(self, Nnut, Epochs=900): #2000
        # self.Nu = Nu
        self.Nnut = Nnut
        self.Epochs = Epochs

    def read(self):
        ## Read the coefficients of Nut
        nut_data_set = np.load("ITHACAoutput/NN/coeffs/coeffL2Nut.npy")
        self.nut_data_set = nut_data_set
        X, Y = self.CreateDataset()
        XTrain, XTest, yTrain,  yTest = train_test_split(
        X, Y, train_size=train_size, random_state=1234, shuffle=False)
        ### Prepare dataset

        scaling = {"scaler_inp": preprocessing.MinMaxScaler( feature_range = (-1, 1) ), "scaler_out": preprocessing.MinMaxScaler(feature_range = (-1, 1) ) } ##feature_range = (-1, 1)

        XTrain = scaling["scaler_inp"].fit_transform(XTrain)

        yTrain = scaling["scaler_out"].fit_transform(yTrain)

        XTest = scaling["scaler_inp"].transform(XTest)

        yTest = scaling["scaler_out"].transform(yTest)

        self.XTrain = XTrain[:, 0:self.Nnut]
        self.yTrain = yTrain[:, 0:self.Nnut]
        self.XTest = XTest[:, 0:self.Nnut]
        self.yTest = yTest[:, 0:self.Nnut]
        self.scaling = scaling

    def CreateDataset(self):
        '''
        : Transform a time series into a prediction dataset
        :  Args: dataset: A numpy ndarray of time series, first dimension is the time steps
        '''
        X, y = [], []
        for i in range(len(self.nut_data_set)-1):
            feature = self.nut_data_set[i]
            target = self.nut_data_set[i+1]
            X.append(feature)
            y.append(target)
        return np.array(X), np.array(y)

    def train_loop(self, dataloader, model, criterion, optimizer):
        #size = len(dataloader.dataset)
        #train_plot = np.ones_like(self.nut_data_set) * np.nan
        # Set the model to training mode - important for batch normalization and dropout layers
        # Unnecessary in this situation but added for best practices

        train_batch_losses = []
        model.train()
        for batch, (X, y) in enumerate(dataloader):
            # Compute prediction and loss
            train_pred = model(X.unsqueeze(0))
            train_loss = criterion(train_pred, y)
            # Backpropagation
            train_loss.backward()
            optimizer.step()
            optimizer.zero_grad()
            #train_loss, current = train_loss.item(), (batch + 1) * len(X)
            train_batch_losses.append(train_loss.detach().numpy())
        loss = np.mean(train_batch_losses)    
        print(f"Train loss: {loss:>10.9f}")
        return loss

    def test_loop(self, testloader, model, criterion):
        # Set the model to evaluation mode - important for batch normalization and dropout layers
        # Unnecessary in this situation but added for best practices

        test_batch_losses = []
        model.eval()
        #size = len(testloader.dataset)
        #num_batches = len(testloader)
        #test_loss, correct = 0.0, 0.0
        #test_plot = np.ones_like(self.nut_data_set) * np.nan
        # Evaluating the model with torch.no_grad() ensures that no gradients are computed during test mode
        # also serves to reduce unnecessary gradient computations and memory usage for tensors with requires_grad=True
        with torch.no_grad():
             for X, y in testloader:
                 pred = model(X.unsqueeze(0))
                 test_loss = criterion(pred, y).item()
                 test_batch_losses.append(test_loss)    
             # correct += (pred.argmax(1) == y).type(torch.float).sum().item()
             #test_loss /= num_batches
             loss  = np.mean(test_batch_losses)    
             print(f"Test loss:  Avg loss: {loss:>10.9f} \n")
             return loss

    def training_testing(self, l_r=1e-4):  #1e-3 #1e-2
        ### Create the model
        model = lstm_seq2seq(self.Nnut, hidden_size)
        #print(model)
        #sys.exit(0)
        self.model = model

        TrainDataset = NutDataset(self.XTrain, self.yTrain)
        self.trainloader = torch.utils.data.DataLoader(
            dataset=TrainDataset, batch_size=batch_size, shuffle=True)
        ### Testing part
        TestDataset = NutDataset(self.XTest, self.yTest)
        self.testloader = torch.utils.data.DataLoader(
            dataset=TestDataset, batch_size=batch_size, shuffle=True)
        ### Optimizer
        optimizer = torch.optim.Adam(model.parameters(), lr=l_r)  # Adam optimizer
        criterion = torch.nn.MSELoss()
        train_plot = []
        test_plot = []
        for t in range(self.Epochs):
            if t%5 ==0:
                print(f"Epoch {t+1}\n-------------------------------")
                train_loss = self.train_loop(
                    self.trainloader, model, criterion, optimizer)
                train_plot.append(train_loss)    
                test_loss = self.test_loop(self.testloader, model, criterion)
                test_plot.append(test_loss)   
        print("Done!")

        plt.plot(train_plot, label="Train_loss")
        plt.plot(test_plot, label="Test_loss")
        plt.legend()
        #plt.autoscale(axis='x',tight=True)
        #plt.autoscale(axis='y',tight=True)
        plt.tight_layout()
        plt.legend(loc='best')
        plt.savefig("Train_test_loss.jpeg")
        plt.show()
        return train_plot, test_plot, model

    def train_test(self):
        self.train_plot, self.test_plot, self.model = self.training_testing()
    def ThePlot(self):

        Test = torch.from_numpy(self.XTest).type(torch.float32)
        test_plot = np.ones_like(self.nut_data_set) * np.nan
        train_plot = np.ones_like(self.nut_data_set) * np.nan

        test_plot[int(len(self.nut_data_set)*train_size) +
           1:len(self.nut_data_set)] = self.model(Test.unsqueeze(0)).detach().numpy()
        XX = torch.from_numpy(self.XTrain).type(torch.float32)
        train_plot[:int(len(self.nut_data_set)*train_size)
              ] = self.model(XX.unsqueeze(0)).detach().numpy() 
              
        #reference = (self.nut_data_set - np.min(self.nut_data_set, axis=0)) /(np.max(self.nut_data_set, axis=0)-np.min(self.nut_data_set, axis=0)) 
        scaler = preprocessing.MinMaxScaler(feature_range = (-1, 1))    
        reference = scaler.fit_transform(self.nut_data_set)       
        ### plot
        f = plt.figure(figsize=(15, 8))
        f.set_tight_layout(True)
        ax1 = f.add_subplot(221)
        ax2 = f.add_subplot(222)
        ax3 = f.add_subplot(223)
        ax4 = f.add_subplot(224)
        #ax5 = f.add_subplot(335)
        #ax6 = f.add_subplot(336)
        #ax7 = f.add_subplot(337)
        #ax8 = f.add_subplot(338)
        #ax9 = f.add_subplot(339)
        #axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]
        axes = [ax1, ax2, ax3, ax4]
        for i in range(4):
            #i = i%4
            axes[i].plot(reference[:, i], 'ko-',
                         label='POD-Gal' + str(i+1), linewidth=2.5)
            axes[i].plot(train_plot[:, i],  'rp--',     label='Training',
                        linewidth=2.5, markerfacecolor='r', alpha=1)
            axes[i].plot(test_plot[:, i],   'bX--',     label='Prediction',
                         linewidth=2.5, markerfacecolor='b', alpha=1)
            axes[i].legend(loc='center')
            #axes[i].autoscale(axis='x',tight=True)
            # axes[i].show()
        #f.legend(bbox_to_anchor =(0.75, 1.15), ncol = 3)
        f.set_tight_layout(True)
        f.savefig("LstmEnDecEddyViscosity.jpeg")   
        plt.show()

    def save(self):
        m = torch.jit.script(self.model)
        print(m.code)
        #sys.exit(0)
        np.save("ITHACAoutput/NN/LstmMinInp_" + str(self.Nnut) + ".npy",self.scaling["scaler_inp"].min_[:,None])
        np.save("ITHACAoutput/NN/LstmScaleInp_"  +  str(self.Nnut) + ".npy",self.scaling["scaler_inp"].scale_[:,None])
        np.save("ITHACAoutput/NN/LstmMinOut_"   + str(self.Nnut) + ".npy",self.scaling["scaler_out"].min_[:,None])
        np.save("ITHACAoutput/NN/LstmScaleOut_"   + str(self.Nnut) + ".npy",self.scaling["scaler_out"].scale_[:,None])
        m.save("ITHACAoutput/NN/LstmNet_" + str(self.Nnut)+".pt")
          
