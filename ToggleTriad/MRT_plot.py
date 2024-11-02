import matplotlib.pyplot as plt
import numpy as np
import re

noises = np.arange(0.0, 0.011,0.002).round(4)
states = {
          '100' : [],
          '010' : [],
          '001' : [],
          '110' : [],
          '011' : [],
          '101' : [],
          '111' : [],
          '000' : []
            }
param = '0'
file = 'stateMRTs_param' + str(param) + '.txt'
f = open("C:/project/Summer_2022_DrMKJolly/Link_Noise/individualSimuls/ToggleTriad/test/" + file, 'r')

for row in f:
  row = re.split(r'\s{1,}', row)
  row = row[:-1]
  if str(row[0]) in list(states.keys()):
    print(row)
    states[row[0]].append(round(float(row[-1]),2))
f.close()


f = plt.figure()
f.set_figwidth(7)
f.set_figheight(5)
for key in states.keys():
    y = states[key]
    x = noises
    plt.plot(x,y, label = key, marker='o', alpha = 0.7)
plt.title('MRTs for TT; Param: ' + str(param))
plt.ylim(-0.1,1500)
plt.yscale('symlog')
#plt.legend(loc = 4)
plt.legend(bbox_to_anchor=(1.04, 1), borderaxespad=0)
#plt.show()
plt.savefig("C:/project/Summer_2022_DrMKJolly/Link_Noise/individualSimuls/ToggleTriad/test/MRTs/" + "MRTs_TT_param_" + str(param) +'.png',bbox_inches="tight", dpi = 300)
plt.close()