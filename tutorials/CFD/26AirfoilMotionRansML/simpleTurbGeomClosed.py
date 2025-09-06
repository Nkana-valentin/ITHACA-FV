import numpy as np
import torch
from torch.autograd import Variable
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.model_selection import train_test_split 
from torch.utils.data import Dataset
from scipy.io import mmread
import math
import sys

torch.manual_seed(0)
torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False
np.random.seed(0)

device = torch.device('cpu')
batch_size = 8
learning_rate =0.00001
weight_decay = 1e-9

plt.style.use('seaborn-whitegrid')
class NutDataset(Dataset):
    def __init__(self, X, y):
        self.X = torch.from_numpy(X).type(torch.float32)
        self.y = torch.from_numpy(y).type(torch.float32)
    def __len__(self):
        return len(self.X)
    def __getitem__(self, idx):
        return [self.X[idx], self.y[idx]]
       

class Net:
    def __init__(self,Nu,Nnut,Epochs):
        self.Nu = Nu
        self.Nnut = Nnut
        self.Epochs = Epochs
        self.batch_size = batch_size
        self.learning_rate = learning_rate
        self.weight_decay = weight_decay
    def read(self):
        ## Read the coefficients train
        # U
        U_data_set = np.load("ITHACAoutput/NN/coeffs/coeffL2U.npy")
        ## Nut
        nut_data_set = np.load("ITHACAoutput/NN/coeffs/coeffL2Nut.npy")
        inp_np_train_U, inp_np_test_U, out_np_train, out_np_test = train_test_split(U_data_set, nut_data_set, test_size=0.3, random_state=1234) 
        
        ## Prepare dataset
        self.inp_np_train = inp_np_train_U[:,0:self.Nu]
        self.inp_np_test =  inp_np_test_U[:,0:self.Nu]
       
        self.out_np_train = out_np_train[:,0:self.Nnut]
        self.out_np_test = out_np_test[:,0:self.Nnut]


    def Net1(self,inp_np_train, inp_np_test, out_np_train, out_np_test, Epochs, learning_rate=learning_rate, weight_decay = weight_decay, batch_size=batch_size):
        # Create NN Network 1

        Nin = inp_np_train.shape[1]
        Nout = out_np_train.shape[1]
        ## Model with forward pass included
        model = torch.nn.Sequential(
             torch.nn.Linear(Nin, 256 ),
             torch.nn.LeakyReLU(),
             torch.nn.Linear(256, 256),
             torch.nn.LeakyReLU(),
             torch.nn.Linear(256, 256),
             torch.nn.LeakyReLU(),
             torch.nn.Linear(256, Nout),
         )
        
        
        
        scaling = {"scaler_inp": preprocessing.MinMaxScaler(feature_range =(-1, 1) ), "scaler_out": preprocessing.MinMaxScaler(feature_range =(-1, 1) )}
        #print("================ inside Net method line 85 =================") 
        inp_np_train = scaling["scaler_inp"].fit_transform(inp_np_train)
        #print("================ inside Net method line 87 =================") 
        out_np_train = scaling["scaler_out"].fit_transform(out_np_train)
        #print("================ inside Net method line 89 =================") 
        inp_np_test = scaling["scaler_inp"].transform(inp_np_test)
        #print("================ inside Net method line 91 =================") 
        out_np_test = scaling["scaler_out"].transform(out_np_test)
        ### Train and test sets of the nut Data set.
        trainset = NutDataset(inp_np_train, out_np_train)
        testset = NutDataset(inp_np_test, out_np_test)
        
        trainloader = torch.utils.data.DataLoader(dataset=trainset, 
                                                  batch_size=batch_size, 
                                                  shuffle=True)
        #trainloader = iter(trainloader)
        testloader = torch.utils.data.DataLoader(dataset=testset, 
                                                 batch_size=batch_size, 
                                                 shuffle=True)
            
        # Loss Functions
        loss_fn = torch.nn.MSELoss()
        optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, weight_decay=weight_decay)

        tplot = []
        lossplottrain = []
        lossplottest = []
        total_samples = len(trainloader)
        print("len(trainloader)", len(trainloader) )
        print("len(testloader)", len(testloader) )
        #number_iterations = total_samples/batch_size
        
        for t in range(self.Epochs):
            # Before the backward pass, use the optimizer object to zero all of the
            # gradients for the variables it will update (which are the learnable
            # weights of the model). This is because by default, gradients are
            # accumulated in buffers( i.e, not overwritten) whenever .backward()
            # is called. Checkout docs of torch.autograd.backward for more details.
            batch_losses = []
            for inputs, labels in trainloader:
                inputs, labels =  inputs.to(device), labels.to(device)
                optimizer.zero_grad() #To empty the gradients 
                # forward + backward + optimize
                outputs = model(inputs) #For prediction
                loss = loss_fn(outputs, labels)
                #loss2 = self.weighted_loss(outputs,labels)
                loss.backward()
                optimizer.step() # to update the weights
                batch_losses.append(loss.item())
            loss = np.mean(batch_losses)
            # evaluate accuracy on test set
            batch_test_losses = []
            model.eval()
            for inputs_test, labels_test in testloader:
                inputs_test, labels_test =  inputs_test.to(device), labels_test.to(device)
                outputs_test = model(inputs_test)
                test_loss = loss_fn(outputs_test, labels_test)
                batch_test_losses.append(test_loss.item())
            test_loss = np.mean(batch_test_losses)
            if t%10 == 9:
                print(t, "loss on train" , loss)
                print(t, "loss on test" , test_loss)
                tplot.append(t)
                lossplottrain.append(loss)
                lossplottest.append(test_loss)
            model.train()

        return tplot,lossplottrain,lossplottest,model,scaling 

    def train(self):
        #print("========== inside train method =========================") 
        self.t_plot, self.lossplottrain, self.lossplottest,self.model,self.scaling = self.Net1(self.inp_np_train, self.inp_np_test, self.out_np_train, self.out_np_test,
                                                                      self.Epochs, self.learning_rate, self.weight_decay, self.batch_size)

                                                                      
    def plot_loss(self):
        plt.plot(self.t_plot, self.lossplottrain, 'ro-',  label="train")
        plt.plot(self.t_plot, self.lossplottest, 'bX-',  label="test")
        plt.tight_layout()
        plt.legend()
        plt.savefig("losses.jpeg")
        plt.show()
    
    def save(self):
        m = torch.jit.script(self.model)
        np.save("ITHACAoutput/NN/minInp_"+str(self.Nu) + "_" +str(self.Nnut) + ".npy",self.scaling["scaler_inp"].min_[:,None])
        np.save("ITHACAoutput/NN/scaleInp_"+str(self.Nu) + "_" +str(self.Nnut) + ".npy",self.scaling["scaler_inp"].scale_[:,None])
        np.save("ITHACAoutput/NN/minOut_"+str(self.Nu) + "_" +str(self.Nnut) + ".npy",self.scaling["scaler_out"].min_[:,None])
        np.save("ITHACAoutput/NN/scaleOut_"+str(self.Nu) + "_" +str(self.Nnut) + ".npy",self.scaling["scaler_out"].scale_[:,None])
        
        m.save("ITHACAoutput/NN/Net_"+str(self.Nu) + "_" +str(self.Nnut)+".pt")
        np.save("ITHACAoutput/NN/trainLoss_"+str(self.Nu) + "_" +str(self.Nnut)+ ".npy",self.lossplottrain)
        np.save("ITHACAoutput/NN/testLoss_"+str(self.Nu) + "_" +str(self.Nnut)+ ".npy",self.lossplottest)
       

