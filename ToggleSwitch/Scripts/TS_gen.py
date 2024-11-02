# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 12:30:48 2022

@author: abhay
"""

import numpy as np
import os
import random
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import statistics as stat

def H(b,b0a,nba):
    return 1 / (1 + (b / b0a) ** nba)

def HS(b, b0a, nba, lba):
    return H(b, b0a, nba) + lba * (1 - H(b, b0a, nba))

def step(x, tau):
    return 1 * (x > tau)

def count_switches(lst):
    switch_count = 0
    for i in range(1, len(lst)):
        if lst[i] != lst[i - 1]:
            switch_count += 1
    return switch_count

def integration_const(p, time, time2, index, p_index, noise):
    dt = time[1] - time[0]
    points = time.size - 1

    st_10 = 0
    st_01 = 0
    st_00 = 0
    st_11 = 0

    curr_st = []
    
    D = np.empty(points + 1)
    E = np.empty(points + 1)

    D[0] = p['Dic'][index]
    E[0] = p['Eic'][index]

    lEtoD = p['lEtoD'][p_index]
    lDtoE = p['lDtoE'][p_index]
    
    gD = p['gD'][p_index]
    gE = p['gE'][p_index]
    
    kD = p['kD'][p_index]
    kE = p['kE'][p_index]
    
    trdEtoD = p['E0'][p_index]
    trdDtoE = p['D0'][p_index]
    
    nEtoD = p['nEtoD'][p_index]
    nDtoE = p['nDtoE'][p_index]
    
    thresh_D = 15.474743669628065#p['D0'][p_index]#(gD / (2*(1+kD))) 
    thresh_E = 17.00017918896313#p['E0'][p_index]#(gE / (2*(1+kE))) 

    lamda_track = [[lEtoD],[lDtoE]]
    for i in range(1, points + 1):
        #if time[i] % 1000 == 0:
        #    print(time[i])

        if time[i] in time2:
            lEtoD = max(0,min(1,(lEtoD + random.gauss(0,noise))))
            lDtoE = max(0,min(1,(lDtoE + random.gauss(0,noise))))

        D[i] = D[i-1] + dt * (gD * HS(E[i - 1], trdEtoD, nEtoD, lEtoD) - kD * D[i-1])
        E[i] = E[i-1] + dt * (gE * HS(D[i - 1], trdDtoE, nDtoE, lDtoE) - kE * E[i-1])
        
        
        if (time[i] > 50.0):
            if D[i] > thresh_D and E[i] > thresh_E:
                st_11 += 0.01
                curr_st.append('st_11')
            elif D[i] > thresh_D and E[i] < thresh_E:
                st_10 += 0.01
                curr_st.append('st_10')
            elif D[i] < thresh_D and E[i] > thresh_E:
                st_01 += 0.01
                curr_st.append('st_01')
            elif D[i] < thresh_D and E[i] < thresh_E:
                st_00 += 0.01
                curr_st.append('st_00')

    num_switches = count_switches(curr_st)

    return D,E,st_10,st_01,st_11,st_00, lamda_track, num_switches, curr_st


T = 1000
dt = 0.01
time = np.arange(0.0,T+dt,dt).round(2)
dt2 = [0.01, 0.1, 1, 10, 100]

noises = [0.006]#np.arange(0.0, 0.011,0.002).round(4)
n_points = time.size

folder = "C:/project/Summer_2022_DrMKJolly/Link_Noise/individualSimuls/ToggleSwitch"

p = {}

#production paramters
p['gD'] = np.array([68.569293, 58.039657, 82.550753, 41.651556, 71.488478, 40.327989, 76.996154, 60.712655, 92.712464, 86.256235, 29.220021, 96.414856, 90.011694, 40.637466, 99.869611, 90.472818, 49.357704, 61.386973, 39.971596, 77.730182, 58.452978, 1.383607, 34.689273, 12.256667, 32.210048, 82.138333, 97.678877])#[91.947511, 50.212433, 64.102184, 83.975461, 39.691404] 86.507167, 86.507167, 53.428293
p['gE'] = np.array([74.210075, 58.18601, 98.508194, 82.175921, 37.186424, 45.131359, 77.485749, 10.343125, 73.861576, 13.277899, 21.498043, 20.956543, 83.591827, 98.165826, 89.769736, 49.995145, 84.77402, 19.619702, 71.123128, 39.907116, 60.956545, 17.90716, 7.93537, 18.154011, 69.224106, 42.201299, 50.140873])#[38.696894, 95.256533, 25.889225, 66.912087, 77.538958] 56.780515, 56.780515, 41.499543

#Threshold constants
p['D0'] = np.array([8.065301, 9.261068, 0.83205, 11.739399, 5.016349, 4.714826, 7.017118, 24.263543, 4.929993, 12.352278, 9.81302, 7.392951, 31.036897, 4.688775, 44.031103, 45.19638, 28.051931, 41.285941, 26.283163, 35.758416, 16.561663, 7.704483, 48.188933, 47.562544, 2.783271, 39.928634, 4.960733])#[53.920926, 14.888436, 31.018099, 53.370193,  3.534087] 35.283372, 35.283372, 34.244496
p['E0'] = np.array([10.718555, 47.958451, 3.84172, 23.868801, 28.27149, 36.045379, 7.576657, 4.191268, 8.331926, 2.676052, 6.754302, 3.05262, 33.466578, 17.891462, 31.479368, 2.219954, 18.265321, 2.843591, 10.792172, 9.547084, 16.351367, 39.843276, 37.67959, 24.646821, 18.280697, 6.593724, 2.613904])#[36.738347, 35.491987,  7.507220, 10.120514, 42.641427] 15.722092, 15.722092, 31.678905

#Degradation constants
p['kD'] = np.array([0.544254, 0.183888, 0.731622, 0.12611, 0.326093, 0.786342, 0.547334, 0.278875, 0.204302, 0.442129, 0.23331, 0.631857, 0.230945, 0.66882, 0.218125, 0.360391, 0.427112, 0.649088, 0.261615, 0.392401, 0.774377, 0.893832, 0.75009, 0.567152, 0.288836, 0.182606, 0.490632])#[0.699470, 0.390561, 0.829410, 0.207175, 0.881824] 0.673448, 0.673448, 0.9119
p['kE'] = np.array([0.257111, 0.320929, 0.867287, 0.388199, 0.565778, 0.541437, 0.822748, 0.222578, 0.211754, 0.856726, 0.503249, 0.99385, 0.417427, 0.65878, 0.343341, 0.285938, 0.207467, 0.300631, 0.750361, 0.519045, 0.250406, 0.950487, 0.883725, 0.885848, 0.334602, 0.649016, 0.830403])#[0.442699, 0.444709, 0.801409, 0.498639, 0.431325] 0.785925, 0.785925, 0.824635

#Fold Change
p['lDtoE'] = np.array([0.061184, 0.296, 0.133122, 0.929518, 0.798589, 0.78079, 0.041854, 0.021802, 0.016178, 0.016149, 0.02563, 0.015496, 0.268771, 0.026455, 0.03954, 0.012613, 0.031571, 0.021592, 0.031398, 0.039624, 0.015844, 0.012172, 0.04193, 0.015208, 0.970356, 0.978338, 0.301517])#[0.011075, 0.028405, 0.017531, 0.015665, 0.011402] 0.068019, 0.068019, 0.034747
p['lEtoD'] = np.array([0.691418, 0.451505, 0.011709, 0.015318, 0.022299, 0.010157, 0.948832, 0.606836, 0.299244, 0.043031, 0.038002, 0.014737, 0.015195, 0.010361, 0.027391, 0.230602, 0.717117, 0.501457, 0.860329, 0.196282, 0.306813, 0.063278, 0.019581, 0.01173, 0.013856, 0.010985, 0.014479])#[0.021791, 0.018186, 0.018946, 0.020940, 0.018222] 0.01201, 0.01201, 0.014156

#Hill Coeff
p['nDtoE'] = np.array([10.0, 8.0, 8.0, 3.0, 1.0, 9.0, 9.0, 7.0, 6.0, 5.0, 8.0, 5.0, 1.0, 1.0, 1.0, 2.0, 9.0, 6.0, 8.0, 1.0, 2.0, 1.0, 10.0, 6.0, 7.0, 9.0, 6.0])#[3.000000, 6.000000, 6.000000, 2.000000, 5.000000] 7.0, 7.0, 3.0
p['nEtoD'] = np.array([10.0, 1.0, 8.0, 8.0, 3.0, 6.0, 2.0, 10.0, 7.0, 9.0, 6.0, 5.0, 2.0, 8.0, 3.0, 10.0, 7.0, 6.0, 6.0, 9.0, 10.0, 9.0, 6.0, 9.0, 8.0, 4.0, 5.0])#[5.000000, 3.000000, 3.000000, 5.000000, 6.000000] 6.0, 6.0, 6.0

params_no_ls = ['0(11)', '2(01)', '3(10)', '4(01_10)', '5(01)', '6(10)', '7(10)', '8(00)', '1(01)']
param_no = '1(01)'

n = 25


for i in range(11, 12, 3):
    param_no = params_no_ls[int(i//3)] + "_log2Snap"
    switching_events = {}

    if os.path.exists(folder + "/Figures/param_no_" + param_no + "/npy_files/size_50_random"+str(i%3)+'.npy'):
        arr = np.load(folder + "/Figures/param_no_" + param_no + "/npy_files/size_50_random"+str(i%3)+'.npy')
        p['Dic'] = arr[0]
        p['Eic'] = arr[1]
    else:
        #Initial Conditions
        p['Dic'] = np.random.random(size = n) * (1.5 * (p['gD'][i]/p['kD'][i]))
        p['Eic'] = np.random.random(size = n) * (1.5 * (p['gE'][i]/p['kE'][i]))

        arr = np.array([p['Dic'],p['Eic']])
        
        if os.path.exists(folder + "/Figures/param_no_" + param_no + "/npy_files") == False:
            os.mkdir(folder + "/Figures/param_no_" + param_no)
            os.mkdir(folder + "/Figures/param_no_" + param_no + "/npy_files")

        np.save(folder + "/Figures/param_no_" + param_no + "/npy_files/size_50_random"+str(int(i%3))+'.npy', arr)
    
    if os.path.exists(folder + "/test/param_no_" + param_no) == False:
        os.mkdir(folder + "/test/param_no_" + param_no)
        
    
    f_mrt = open(folder + '/test/param_no_' + param_no + '/stateMRTs_param' + str(i) + '.txt','w')
    for timestep in dt2:
        time2 = np.arange(50.0,T+dt,timestep).round(2)
        f_mrt.write('\nFor timestep: {} and parameter: {}'.format(timestep, i))
        for noise in noises:
            tt_dynamics = []
            
            states = {
                '10' : [],
                '01' : [],
                '00' : [],
                '11' : []
                }

            total_switches_ls = []
            f = plt.figure(figsize = (17, 15))
            ax1 = plt.subplot(331)
            ax2 = plt.subplot(332)
            ax3 = plt.subplot(333)
            ax4 = plt.subplot(312)
            ax5 = plt.subplot(313)
            ax = [ax1,ax2,ax3,ax4,ax5]
            
            for l in range(n):
                D,E,st_10,st_01,st_11,st_00,lt, num_sw, curr_st = integration_const(p, time, time2, l, i, noise)
                total_switches_ls.append(num_sw)
                
                states['10'].append(st_10)
                states['01'].append(st_01)
                states['00'].append(st_00)
                states['11'].append(st_11)
                
                ax[0].plot(time, D, linewidth = 1.5, label = str(p['Dic'][l]))
                ax[0].set_ylabel("A",fontsize=14)
                ax[1].plot(time, E, linewidth = 1.5, label = str(p['Eic'][l]))
                ax[1].set_ylabel("B",fontsize=14)
                
                ax[2].plot(time, D, linewidth = 1.5, label = 'A', color = 'red', alpha = 0.1)
                ax[2].plot(time, E, linewidth = 1.5, label = 'B', color = 'blue', alpha = 0.1)
                ax[2].set_ylabel("AB",fontsize=14)

                ax[3].plot(time[5001:], curr_st)
                ax[4].plot(time[5001 + 20000:5001 + 30001], curr_st[20000:30001])
                
                tt_dynamics.append(np.array([D,E]))
            
            switching_events[noise] = [sum(total_switches_ls),stat.stdev(total_switches_ls),stat.mean(total_switches_ls)]
            
            
            f_mrt.write('\nFor noise: {}'.format(noise))
            f_mrt.write('\n#SE: {:<8} {:<15} {:<10}'.format(switching_events[noise][0], switching_events[noise][1], switching_events[noise][2]))
            f_mrt.write('\n{:<8} {:<15}'.format('State', 'MRT'))
            f_mrt.write('\n{:<8} {:<15} {:<10}'.format('10', stat.stdev(states['10']), stat.mean(states['10'])))
            f_mrt.write('\n{:<8} {:<15} {:<10}'.format('01', stat.stdev(states['01']), stat.mean(states['01'])))
            f_mrt.write('\n{:<8} {:<15} {:<10}'.format('00', stat.stdev(states['00']), stat.mean(states['00'])))
            f_mrt.write('\n{:<8} {:<15} {:<10}'.format('11', stat.stdev(states['11']), stat.mean(states['11'])))

            
            
            plt.tick_params(axis='y',labelsize=12,rotation=90)
            plt.tick_params(axis='x',labelsize=12)
            f.suptitle('TS_Control')
            ax[1].set_xlabel("Time",fontsize=14)
            
            red_patch = mpatches.Patch(color='red', label='A')
            blue_patch = mpatches.Patch(color='blue', label='B')
            
            ax[2].legend(handles = [red_patch, blue_patch], loc = 'center right')
            
            f.tight_layout()
            #plt.show()
            plt.savefig(folder + "/test/param_no_" + param_no + "/Snaptime_dynamics_Noise" + str(noise) + "_dt2" + str(timestep) + "_param_" + str(int(i%3)) +'.png',dpi = 300)
            plt.close()        
            #arr_dynamics.append(tt_dynamics)

            
                
            #file = 'value_storage/toggleTriad_noise_' + str(noise) + '_time_' + str(T) + '_param_' + str(i) + '.npy'
            #np.save(folder + file, np.array(tt_dynamics))


    '''
        f = plt.figure()
        f.set_figwidth(7)
        f.set_figheight(5)
        
        plt.bar(range(len(switching_events)), [i[0] for i in list(switching_events.values())], align='center')
        #plt.errorbar(range(len(switching_events)), [i[2] for i in list(switching_events.values())], yerr = [i[1] for i in list(switching_events.values())])
        plt.xticks(range(len(switching_events)), list(switching_events.keys()))
        

        plt.ylabel('# Switching events')
        plt.xlabel('Noise')
        plt.title("Num switching evenets; Param: " + str(i))

        plt.savefig(folder + "/test/param_no_" + param_no + "/TS_Num_switching_param_" + str(int(i%3)) +'.png', dpi = 300)
        plt.close()

        f = plt.figure()
        f.set_figwidth(7)
        f.set_figheight(5)
        
        plt.bar(range(len(switching_events)), [i[2] for i in list(switching_events.values())], align='center')
        plt.errorbar(range(len(switching_events)), [i[2] for i in list(switching_events.values())], yerr = [i[1] for i in list(switching_events.values())], fmt = 'o', color = 'r')
        plt.xticks(range(len(switching_events)), list(switching_events.keys()))
        

        plt.ylabel('# Switching events')
        plt.xlabel('Noise')
        plt.title("Num switching evenets; Param: " + str(i))
        plt.savefig(folder + "/test/param_no_" + param_no + "/TS_Num_switching_error_param_" + str(int(i%3)) +'.png',dpi = 300)

        #plt.savefig("C:/project/Summer_2022_DrMKJolly/Link_Noise/individualSimuls/ToggleTriad/MRTs/" + "test_TT_Num_switching_error_param_" + str(i) +'.png')
        plt.close()
    '''
    f_mrt.close()
    print('file '+str(i)+'saved')