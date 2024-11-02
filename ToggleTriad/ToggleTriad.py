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
    
    st_100 = 0
    st_010 = 0
    st_001 = 0
    st_110 = 0
    st_011 = 0
    st_101 = 0
    st_111 = 0    
    st_000 = 0

    curr_state = []

    A = np.empty(points + 1)
    B = np.empty(points + 1)
    C = np.empty(points + 1)
    
    A[0] = p['Aic'][index]
    B[0] = p['Bic'][index]
    C[0] = p['Cic'][index]

    lAtoB = p['lAtoB'][p_index]
    lBtoC = p['lBtoC'][p_index]
    lCtoA = p['lCtoA'][p_index]
    lAtoC = p['lAtoC'][p_index]
    lCtoB = p['lCtoB'][p_index]
    lBtoA = p['lBtoA'][p_index]
    
    ga = p['gA'][p_index]
    gb = p['gB'][p_index]
    gc = p['gC'][p_index]
    
    ka = p['kA'][p_index]
    kb = p['kB'][p_index]
    kc = p['kC'][p_index]
    
    trdCtoA = p['trdCtoA'][p_index]
    trdBtoA = p['trdBtoA'][p_index]
    trdAtoB = p['trdAtoB'][p_index]
    trdCtoB = p['trdCtoB'][p_index]
    trdAtoC = p['trdAtoC'][p_index]
    trdBtoC = p['trdBtoC'][p_index]
    
    nCtoA = p['nCtoA'][p_index]
    nBtoA = p['nBtoA'][p_index]
    nAtoB = p['nAtoB'][p_index]
    nCtoB = p['nCtoB'][p_index]
    nAtoC = p['nAtoC'][p_index]
    nBtoC = p['nBtoC'][p_index]

    threshA = 4.922509060330396#(ga / ka) / 2
    threshB = 4.954690196044068#(gb / kb) / 2
    threshC = 5.308042259014317#(gc / kc) / 2

    for i in range(1, points + 1):
        #if time[i] % 1000 == 0:
        #    print(time[i])
        
        if time[i] in time2:
            lAtoB = max(0,min(1,(lAtoB + random.gauss(0,noise))))
            lBtoC = max(0,min(1,(lBtoC + random.gauss(0,noise))))
            lCtoA = max(0,min(1,(lCtoA + random.gauss(0,noise))))
            lAtoC = max(0,min(1,(lAtoC + random.gauss(0,noise))))
            lCtoB = max(0,min(1,(lCtoB + random.gauss(0,noise))))
            lBtoA = max(0,min(1,(lBtoA + random.gauss(0,noise))))
        

        
        A[i] = A[i-1] + dt * (ga * HS(C[i - 1], trdCtoA, nCtoA, lCtoA) * HS(B[i - 1], trdBtoA, nBtoA, lBtoA) - ka * A[i-1])
        B[i] = B[i-1] + dt * (gb * HS(A[i - 1], trdAtoB, nAtoB, lAtoB) * HS(C[i - 1], trdCtoB, nCtoB, lCtoB) - kb * B[i-1])
        C[i] = C[i-1] + dt * (gc * HS(B[i - 1], trdBtoC, nBtoC, lBtoC) * HS(A[i - 1], trdAtoC, nAtoC, lAtoC) - kc * C[i-1])
        

        if (time[i] > 55.0):
            if A[i] > threshA and B[i] < threshB and C[i] < threshC:
                st_100 += 0.01
                curr_state.append('st_100')
            elif A[i] < threshA and B[i] > threshB and C[i] < threshC:
                st_010 += 0.01
                curr_state.append('st_010')
            elif A[i] < threshA and B[i] < threshB and C[i] > threshC:
                st_001 += 0.01
                curr_state.append('st_001')
            elif A[i] > threshA and B[i] > threshB and C[i] < threshC:
                st_110 += 0.01
                curr_state.append('st_110')
            elif A[i] < threshA and B[i] > threshB and C[i] > threshC:
                st_011 += 0.01
                curr_state.append('st_011')
            elif A[i] > threshA and B[i] < threshB and C[i] > threshC:
                st_101 += 0.01
                curr_state.append('st_101')
            elif A[i] > threshA and B[i] > threshB and C[i] > threshC:
                st_111 += 0.01
                curr_state.append('st_111')
            elif A[i] < threshA and B[i] < threshB and C[i] < threshC:
                st_000 += 0.01
                curr_state.append('st_000')

    num_switches = count_switches(curr_state)

    return A,B,C,st_100,st_010,st_001,st_110,st_011,st_101,st_111,st_000, num_switches



T = 1000
dt = 0.01
time = np.arange(0.0,T+dt,dt).round(2)
dt2 = 0.01
time2 = np.arange(50.0,1000+dt,dt2).round(2)
n_points = time.size

noises = np.arange(0.0, 0.011,0.002).round(4)#np.arange(0.0, 0.00051, 0.0001).round(4)

folder = "C:/project/Summer_2022_DrMKJolly/Link_Noise/individualSimuls/ToggleTriad/test/"

p = {}

#production paramters
p['gA'] = [14.311011, 41.559917, 43.215851, 37.135048, 45.60133, 84.229054]#[12.018249, 5.230785, 92.301214, 59.528363, 43.125467, 74.661273]
p['gB'] = [36.830513, 24.881922, 85.037079, 44.572465, 87.445041, 55.375012]#[2.675641, 76.596182, 47.940154, 79.558099, 42.705912, 88.839696]
p['gC'] = [5.767979, 73.833359, 75.617525, 83.198912, 32.735875, 56.632195]#[60.869333, 6.714682, 46.216528, 89.527126, 22.316457, 33.342705]

#Threshold constants
p['trdAtoB'] = [11.549171, 4.776206, 1.286003, 0.212738, 7.601168, 7.853812]#[8.091984 , 3.091889 , 0.400201, 2.367581, 7.893864, 5.979105 ]
p['trdBtoC'] = [5.167491, 5.92531, 12.37715, 13.931096, 9.81062, 13.916002]#[2.561299 , 7.347394 ,  8.68131, 5.630365, 9.875833, 11.509482]
p['trdCtoA'] = [10.563554, 12.741537, 4.079019, 2.296203, 7.040644, 8.988097]#[0.329027 , 6.569456 , 5.929507, 7.94269 , 9.151694, 4.745698 ]
p['trdAtoC'] = [7.745853, 5.466528, 8.059633, 3.40674, 0.543483, 12.10012]#[9.705534 , 9.334428 , 5.253287, 4.583083, 6.194162, 9.969809 ]
p['trdCtoB'] = [4.539267, 8.603515, 0.18813, 2.702384, 2.247522, 0.184975]#[4.111351 , 6.053961 , 7.564845, 3.941795, 12.18131, 8.414342 ]
p['trdBtoA'] = [0.39581, 1.650777, 4.985056, 1.107473, 2.33763, 0.846264]#[11.313256, 10.165773, 6.924967, 7.630631, 9.723269, 11.243532]

#Degradation constants
p['kA'] = [0.416011, 0.727744, 0.741633, 0.881418, 0.290622, 0.842896]#[0.804511, 0.577791, 0.361671, 0.940077, 0.631683, 0.7806  ]
p['kB'] = [0.663866, 0.440636, 0.392892, 0.715639, 0.420955, 0.17808]#[0.178982, 0.534397, 0.949621, 0.509631, 0.758171, 0.913511]
p['kC'] = [0.14316, 0.345527, 0.465843, 0.779907, 0.441988, 0.137671]#[0.767847, 0.180412, 0.230333, 0.899618, 0.429697, 0.678669]

#Fold Change
p['lAtoB'] = [0.013573, 0.055652, 0.013308, 0.011081, 0.012494, 0.020097]#[0.080769, 0.013618, 0.117883, 0.012118, 0.013122, 0.020298]
p['lBtoC'] = [0.030464, 0.012486, 0.026228, 0.012703, 0.024463, 0.014486]#[0.01015 , 0.030059, 0.019566, 0.014407, 0.013806, 0.022222]
p['lCtoA'] = [0.060351, 0.136092, 0.054288, 0.014439, 0.010476, 0.208629]#[0.013822, 0.010811, 0.014476, 0.02025 , 0.017554, 0.012852]
p['lAtoC'] = [0.015798, 0.018996, 0.010448, 0.011854, 0.161245, 0.011876]#[0.013274, 0.021265, 0.012013, 0.010072, 0.011343, 0.010311]
p['lCtoB'] = [0.010536, 0.018326, 0.010579, 0.105287, 0.045982, 0.012697]#[0.017927, 0.031853, 0.011715, 0.011788, 0.069397, 0.016738]
p['lBtoA'] = [0.029469, 0.012697, 0.174123, 0.117494, 0.032876, 0.019761]#[0.017287, 0.05123 , 0.022806, 0.012885, 0.024531, 0.056011]

#Hill Coeff
p['nAtoB'] = [2.0, 5.0, 2.0, 4.0, 1.0, 5.0]#[6.000000, 4.000000, 1.000000, 6.000000, 4.000000, 6.000000]
p['nBtoC'] = [6.0, 6.0, 1.0, 3.0, 6.0, 2.0]#[4.000000, 2.000000, 4.000000, 3.000000, 4.000000, 4.000000]
p['nCtoA'] = [1.0, 2.0, 3.0, 5.0, 5.0, 5.0]#[1.000000, 5.000000, 4.000000, 6.000000, 5.000000, 5.000000]
p['nAtoC'] = [6.0, 1.0, 5.0, 6.0, 3.0, 4.0]#[4.000000, 1.000000, 6.000000, 3.000000, 3.000000, 5.000000]
p['nCtoB'] = [4.0, 4.0, 1.0, 3.0, 6.0, 1.0]#[5.000000, 4.000000, 3.000000, 2.000000, 6.000000, 2.000000]
p['nBtoA'] = [2.0, 2.0, 2.0, 1.0, 1.0, 5.0]#[4.000000, 3.000000, 4.000000, 3.000000, 6.000000, 6.000000]

n = 25

for i in range(0, 6, 1):
    switching_events = {}

    if os.path.exists(folder + 'value_store/size_' + str(n) + '_random' + str(i) + '.npy'):
        arr = np.load(folder + 'value_store/size_' + str(n) + '_random' + str(i) + '.npy')
        p['Aic'] = arr[0]
        p['Bic'] = arr[1]
        p['Cic'] = arr[2]
    else:
        #Initial Conditions
        p['Aic'] = np.random.random(size = n) * (1.5 * (p['gA'][i]/p['kA'][i]))
        p['Bic'] = np.random.random(size = n) * (1.5 * (p['gB'][i]/p['kB'][i]))
        p['Cic'] = np.random.random(size = n) * (1.5 * (p['gC'][i]/p['kC'][i]))

        arr = np.array([p['Aic'],p['Bic'],p['Cic']])

        np.save(folder + 'value_store/size_' + str(n) + '_random' + str(i) + '.npy', arr)
        
    f_mrt = open(folder + 'stateMRTs_param' + str(i) + '.txt','w')
    
    for noise in noises:
        tt_dynamics = []
        
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

        total_switches_ls = []
        f, ax = plt.subplots(1,4,figsize = (20, 5), sharex = True, sharey = False)
        
        for l in range(n):
            A,B,C,st_100,st_010,st_001,st_110,st_011,st_101,st_111,st_000,num_sw = integration_const(p, time, time2, l, i, noise)
            total_switches_ls.append(num_sw)
            
            states['100'].append(st_100)
            states['010'].append(st_010)
            states['001'].append(st_001)
            states['110'].append(st_110)
            states['011'].append(st_011)
            states['101'].append(st_101)
            states['111'].append(st_111)            
            states['000'].append(st_000)
            
            ax[0].plot(time[500:], A[500:], linewidth = 1.5, label = str(p['Aic'][l]))
            ax[0].set_ylabel("A",fontsize=14)
            ax[1].plot(time[500:], B[500:], linewidth = 1.5, label = str(p['Bic'][l]))
            ax[1].set_ylabel("B",fontsize=14)
            ax[2].plot(time[500:], C[500:], linewidth = 1.5, label = str(p['Cic'][l]))
            ax[2].set_ylabel("C",fontsize=14)
            
            ax[3].plot(time[500:], A[500:], linewidth = 1.5, label = 'A', color = 'red', alpha = 0.1)
            ax[3].plot(time[500:], B[500:], linewidth = 1.5, label = 'B', color = 'blue', alpha = 0.1)
            ax[3].plot(time[500:], C[500:], linewidth = 1.5, label = 'C', color = 'green', alpha = 0.1)
            ax[3].set_ylabel("ABC",fontsize=14)
            
            tt_dynamics.append(np.array([A,B,C]))
        
        switching_events[noise] = [sum(total_switches_ls),stat.stdev(total_switches_ls),stat.mean(total_switches_ls)]
        
        
        f_mrt.write('\nFor noise: {} and parameter: {}'.format(noise, i))
        f_mrt.write('\n{:<8} {:<15}'.format('State', 'MRT'))
        f_mrt.write('\n{:<8} {:<15} {:<10}'.format('100', stat.stdev(states['100']), stat.mean(states['100'])))
        f_mrt.write('\n{:<8} {:<15} {:<10}'.format('010', stat.stdev(states['010']), stat.mean(states['010'])))
        f_mrt.write('\n{:<8} {:<15} {:<10}'.format('001', stat.stdev(states['001']), stat.mean(states['001'])))
        f_mrt.write('\n{:<8} {:<15} {:<10}'.format('110', stat.stdev(states['110']), stat.mean(states['110'])))
        f_mrt.write('\n{:<8} {:<15} {:<10}'.format('011', stat.stdev(states['011']), stat.mean(states['011'])))
        f_mrt.write('\n{:<8} {:<15} {:<10}'.format('101', stat.stdev(states['101']), stat.mean(states['101'])))
        f_mrt.write('\n{:<8} {:<15} {:<10}'.format('111', stat.stdev(states['111']), stat.mean(states['111'])))
        f_mrt.write('\n{:<8} {:<15} {:<10}'.format('000', stat.stdev(states['000']), stat.mean(states['000'])))
        
        
        plt.tick_params(axis='y',labelsize=12,rotation=90)
        plt.tick_params(axis='x',labelsize=12)
        f.suptitle('TT_Control')
        ax[1].set_xlabel("Time (in weeks)",fontsize=14)
        
        red_patch = mpatches.Patch(color='red', label='A')
        blue_patch = mpatches.Patch(color='blue', label='B')
        green_patch = mpatches.Patch(color='green', label='C')
        
        ax[3].legend(handles = [red_patch, blue_patch, green_patch], loc = 'center right')
        
        f.tight_layout()
        #plt.show()
        plt.savefig(folder + "ToggleTriad_param_" + str(i) + "_noise_" + str(noise) + ".png",dpi=300)
        plt.close()        
        #arr_dynamics.append(tt_dynamics)

        
            
        #file = 'value_storage/toggleTriad_noise_' + str(noise) + '_time_' + str(T) + '_param_' + str(i) + '.npy'
        #np.save(folder + file, np.array(tt_dynamics))

        
    f_mrt.close()
    print('file '+str(i)+'saved')


    f = plt.figure()
    f.set_figwidth(7)
    f.set_figheight(5)
    
    plt.bar(range(len(switching_events)), [i[0] for i in list(switching_events.values())], align='center')
    #plt.errorbar(range(len(switching_events)), [i[2] for i in list(switching_events.values())], yerr = [i[1] for i in list(switching_events.values())])
    plt.xticks(range(len(switching_events)), list(switching_events.keys()))
    

    plt.ylabel('# Switching events')
    plt.xlabel('Noise')
    plt.title("Num switching evenets; Param: " + str(i))

    plt.savefig(folder + "TT_Num_switching_param_" + str(i) +'.png')
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

    plt.savefig(folder + "TT_Num_switching_error_param_" + str(i) +'.png')
    plt.close()