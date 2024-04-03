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
batch_size = 10 #32
learning_rate = 1e-3  #1e-4 #0.0025 #0.00001 #0.02
weight_decay = 1e-9

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "fantasy",
     "font.size": 24, 
    "font.fantasy": 'Times New Roman',
    'figure.figsize' : (15,5),
    'figure.dpi' : 100
})
plt.style.use('seaborn-whitegrid')
       

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
        ## U
        U_data_set = np.load("ITHACAoutput/NN/coeffs/coeffL2U.npy")
        U_data_set = U_data_set[0:1001]
        self.U_data_set = U_data_set
       
        ## Nut
        nut_data_set = np.load("ITHACAoutput/NN/coeffs/coeffL2Nut.npy")
        nut_data_set = nut_data_set[0:1001]
        self.nut_data_set = nut_data_set
        
        
        inp_np_train_U, inp_np_test_U, out_np_train, out_np_test = train_test_split(U_data_set, nut_data_set, test_size=0.2, random_state=1234,shuffle=False) 
        ## Prepare dataset
        self.inp_np_train = inp_np_train_U[:,0:self.Nu]
        self.inp_np_test =  inp_np_test_U[:,0:self.Nu]
        
        self.out_np_train = out_np_train[:,0:self.Nnut]
        self.out_np_test = out_np_test[:,0:self.Nnut]


    def Net1(self,inp_np_train, inp_np_test, out_np_train, out_np_test, Epochs, learning_rate=learning_rate, batch_size=batch_size):
        # Create NN Network 1

        Nin = inp_np_train.shape[1]
        Nout = out_np_train.shape[1]
        ## Model with forward pass included
        
        model = torch.nn.Sequential(
             torch.nn.Linear(Nin, 256 ),
             torch.nn.ELU(),
             torch.nn.Linear(256, 128),
             torch.nn.ELU(),
             torch.nn.Linear(128, Nout)
         )
        
        scaling = {"scaler_inp": preprocessing.MinMaxScaler(feature_range =(-1, 1) ), "scaler_out": preprocessing.MinMaxScaler(feature_range =(-1, 1) )}
        inp_np_train = scaling["scaler_inp"].fit_transform(inp_np_train)
        out_np_train = scaling["scaler_out"].fit_transform(out_np_train)
        inp_np_test = scaling["scaler_inp"].transform(inp_np_test)
        out_np_test = scaling["scaler_out"].transform(out_np_test)

        
        X_train = torch.tensor(inp_np_train, dtype=torch.float32)
        #print(X_train[0])
        y_train = torch.tensor(out_np_train, dtype=torch.float32)
      
        X_test = torch.tensor(inp_np_test, dtype=torch.float32)
        y_test = torch.tensor(out_np_test, dtype=torch.float32)

            
        # Loss Functions
        loss_fn = torch.nn.MSELoss()
        optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
        #optimizer = torch.optim.SGD(model.parameters(), lr=learning_rate)
        tplot = []
        lossplottrain = []
        lossplottest = []
        
        for epoch in range(Epochs):
            # Before the backward pass, use the optimizer object to zero all of the
            # gradients for the variables it will update (which are the learnable
            # weights of the model). This is because by default, gradients are
            # accumulated in buffers( i.e, not overwritten) whenever .backward()
            # is called. Checkout docs of torch.autograd.backward for more details.
            batch_losses = []
            model.train()
            for i in range(0, len(X_train), batch_size): ####for batch in trainloader:
                Xbatch = X_train[i:i+batch_size]
                y_pred = model(Xbatch)
                ybatch = y_train[i:i+batch_size]      #X, y  = batch
                train_loss = loss_fn(ybatch, y_pred)
                optimizer.zero_grad() 
                train_loss.backward()
                optimizer.step() # to update the weights
            print(f'Finished epoch {epoch}, latest loss {train_loss}')    
                #batch_losses.append(train_loss.item())
            #train_loss = np.mean(batch_losses)
            # evaluate accuracy on test set
            # batch_test_losses = []
            # model.eval()
            # for inputs_test, labels_test in testloader:
            #     #inputs_test, labels_test =  inputs_test.to(device), labels_test.to(device)
            #     outputs_test = model(inputs_test)
            #     #print("outputs_test", outputs_test)
            #     test_loss = loss_fn(outputs_test, labels_test)
            #     batch_test_losses.append(test_loss.item())
            # test_loss = np.mean(batch_test_losses)
            # if epoch%100 == 0:
            #     print(epoch, "train_loss \t" , train_loss) #### "\t", "test_loss \t", test_loss)
            #     tplot.append(epoch)
            #     lossplottrain.append(train_loss)
            #     #lossplottest.append(test_loss)
        ### Compute accuracy
        with torch.no_grad():
            y = model(X_train) ## Validation
            y_pred = model(X_test) ## Validation
            yy = y_pred.detach().numpy()
            zz = y.detach().numpy()
            k = 1
            plt.plot(np.arange(len(zz[:, k]), 1001), yy[:, k], 'b.',linestyle='dashed',  label='$Validation$')
            plt.plot(np.arange(len(zz[:, k])), zz[:, k], 'r.', linestyle='dashed', label='$Truth $'+str(k))
            plt.tight_layout()
            plt.legend(ncol=2, loc='upper left', fontsize=20, bbox_to_anchor=(0.95, -0.25))
            plt.savefig("NNs_" + str(self.Nu) + "_" + str(self.Nnut) +  ".pdf")
            plt.show()
            # print("y_pred = ", y_pred)
            # print("y_train", y_train)
            # acc = y_pred-y_train
            # print(acc)

        return tplot,lossplottrain,lossplottest,model,scaling 

    def train(self):
        #print("========== inside train method =========================") 
        self.t_plot, self.lossplottrain, self.lossplottest,self.model,self.scaling = self.Net1(self.inp_np_train, self.inp_np_test, self.out_np_train, self.out_np_test,
                                                                      self.Epochs, self.learning_rate, self.batch_size)
    def plot_loss(self):
        plt.plot(self.t_plot, self.lossplottrain, 'r--',  label="$Train$", marker='h')
        plt.plot(self.t_plot, self.lossplottest, 'b--',   label="$Validation$",marker='D')
        plt.xlabel("$Epochs$")
        plt.ylabel("$Losses$")
        plt.yscale('log')
        plt.tight_layout()
        plt.legend(ncol=2)
        plt.savefig("losses_" + str(self.Nu) + "_" + str(self.Nnut) +  ".pdf")
        #plt.show()
    
    def save(self):
        m = torch.jit.script(self.model)
        np.save("ITHACAoutput/NN/minInp_"+str(self.Nu) + "_" +str(self.Nnut) + ".npy",self.scaling["scaler_inp"].min_[:,None])
        np.save("ITHACAoutput/NN/scaleInp_"+str(self.Nu) + "_" +str(self.Nnut) + ".npy",self.scaling["scaler_inp"].scale_[:,None])
        np.save("ITHACAoutput/NN/minOut_"+str(self.Nu) + "_" +str(self.Nnut) + ".npy",self.scaling["scaler_out"].min_[:,None])
        np.save("ITHACAoutput/NN/scaleOut_"+str(self.Nu) + "_" +str(self.Nnut) + ".npy",self.scaling["scaler_out"].scale_[:,None])
        
        m.save("ITHACAoutput/NN/Net_"+str(self.Nu) + "_" +str(self.Nnut)+".pt")
        np.save("ITHACAoutput/NN/trainLoss_"+str(self.Nu) + "_" +str(self.Nnut)+ ".npy",self.lossplottrain)
        np.save("ITHACAoutput/NN/testLoss_"+str(self.Nu) + "_" +str(self.Nnut)+ ".npy",self.lossplottest)
       

