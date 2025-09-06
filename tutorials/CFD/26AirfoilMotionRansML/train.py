#from withNet import *
from simpleTurbGeomClosed import *
from files import *

nNut = read_variable("NmodesNutProj", "system/ITHACAdict")
NU = read_variable("NmodesUproj", "system/ITHACAdict")
#NP = read_variable("NmodesPproj", "system/ITHACAdict")

epochs = read_variable("epochs", "system/ITHACAdict")

Netok = Net(NU, nNut, epochs)
Netok.read()

Netok.train()
Netok.plot_loss()
Netok.save()

