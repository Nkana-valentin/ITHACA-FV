#from withNet import *

#from SimpleLstmViscosity import *
from lstmEncoderDecoder import *
from files import *

import sys
#nNut = 15
nNut = read_variable("NmodesNutProj", "system/ITHACAdict")

#epochs = read_variable("epochs", "system/ITHACAdict")

Netok = Net(nNut)
#print(Netok.decoder)
Netok.read()
#sys.exit(0)
Netok.training_testing()
#sys.exit()
Netok.ThePlot()

Netok.save()
