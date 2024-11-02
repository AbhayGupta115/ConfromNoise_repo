import matplotlib.pyplot as plt
import numpy as np
import re

noises = [0.006]#np.arange(0.0, 0.011,0.002).round(4)
dt2 = ['0.01','0.1', '1', '10', '100']
states = {
          '10' : [],
          '01' : [],
          '00' : [],
          '11' : []
            }
switching_events = {}
param = '11'
params_no_ls = ['0(11)', '2(01)', '3(10)', '4(01_10)', '5(01)', '6(10)', '7(10)', '8(00)', '1(01)']
param_no = params_no_ls[int(int(param)//3)]  + "_log2Snap"
file = '/stateMRTs_param' + str(int(param)) + '.txt'
f = open('C:/project/Summer_2022_DrMKJolly/Link_Noise/individualSimuls/ToggleSwitch/test/param_no_' + str(param_no) + file, 'r')

for row in f:
  row = re.split(r'\s{1,}', row)
  row = row[:-1]
  print(row)

  if (len(row) == 6):
    timestep = row[2]
  if (len(row) == 3) and (row[0] == 'For'):
    noise_lvl = row[2]

  if row[0] == '#SE:':
    switching_events[str(noise_lvl) + str(timestep)] = [row[1], row[2], row[3]]

  if str(row[0]) in list(states.keys()):
    
    states[row[0]].append(round(float(row[-1]),2))
f.close()

for i in range(len(dt2)):
  f = plt.figure()
  f.set_figwidth(7)
  f.set_figheight(5)
  for key in states.keys():
      y = states[key][i*len(noises):(i+1)*len(noises)]
      x = noises
      plt.plot(x,y, label = key, marker='o', alpha = 0.7)
  plt.title('MRTs for TS; Param: ' + str(param) + 'timestep:' + str(dt2[i]))
  plt.ylim(-0.1,1500)
  plt.yscale('symlog')
  #plt.legend(loc = 4)
  plt.legend(bbox_to_anchor=(1.15, 1), borderaxespad=0)
  #plt.show()
  plt.savefig('C:/project/Summer_2022_DrMKJolly/Link_Noise/individualSimuls/ToggleSwitch/test/param_no_' + str(param_no) + "/MRTs_TS_param_" + str(int(param)%3) + "_timestep_"+ str(dt2[i]) +'.png',bbox_inches="tight", dpi = 300)
  plt.close()

for i in range(len(noises)):
  f = plt.figure()
  f.set_figwidth(7)
  f.set_figheight(5)
  for key in states.keys():
      y = states[key][i::len(noises)]
      x = dt2
      plt.plot(x,y, label = key, marker='o', alpha = 0.7)
  plt.title('MRTs for TS; Param: ' + str(param) + 'noise:' + str(noises[i]))
  plt.ylim(-0.1,1500)
  plt.yscale('symlog')
  #plt.legend(loc = 4)
  plt.legend(bbox_to_anchor=(1.15, 1), borderaxespad=0)
  #plt.show()
  plt.savefig('C:/project/Summer_2022_DrMKJolly/Link_Noise/individualSimuls/ToggleSwitch/test/param_no_' + str(param_no) + "/MRTs_TS_param_" + str(int(param)%3) + "_noise_"+ str(noises[i]) +'.png',bbox_inches="tight", dpi = 300)
  plt.close()