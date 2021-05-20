import sys
import os
import glob
import numpy as np
from pyMCDS_cells import pyMCDS_cells
import matplotlib.pyplot as plt

argc = len(sys.argv)-1
print("# args=",argc)

#data_dir = 'output'
if (argc < 1):
#  data_dir = int(sys.argv[kdx])
  print("Usage: provide output subdir")
  sys.exit(-1)

kdx = 1
data_dir = sys.argv[kdx]
print('data_dir = ',data_dir)
os.chdir(data_dir)
xml_files = glob.glob('output*.xml')
os.chdir('..')
xml_files.sort()
#print('xml_files = ',xml_files)

ds_count = len(xml_files)
print("----- ds_count = ",ds_count)
mcds = [pyMCDS_cells(xml_files[i], data_dir) for i in range(ds_count)]

tval = np.linspace(0, mcds[-1].get_time(), ds_count)
print('tval= ',tval)

y_energy = np.array( [mcds[idx].data['discrete_cells']['energy'] for idx in range(ds_count)] )
print(y_energy)

plt.plot(tval, y_energy, label='energy', linewidth=1)  #, color='lime')

#plt.legend(loc='center left', prop={'size': 10})

plt.title(data_dir + ": energy")
#plt.savefig(data_dir + '.png')
plt.show()
