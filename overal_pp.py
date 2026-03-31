#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  9 12:57:06 2025

@author: hossein
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
import ast
import matplotlib.animation as animation
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import linregress
import json
import subprocess
from scipy.interpolate import UnivariateSpline
import time as time_package

def beta_a_norm_plotter_C(beta_a_norm_dict, celltype): # Works only for type C
    
    def beta_norm_analysis(x, y):
        
        max_ind = np.argmax(y)
        i_dum = max_ind
        j_dum = max_ind
        
        y_right = y[max_ind:].copy()
        
        integral_goal = np.log(2.0)
        integral_val = 0.0
        
        while integral_val < integral_goal:
            i_dum -= 1
            j_dum = np.argmin(np.abs(y_right-y[i_dum])) + max_ind
            integral_val = dt*np.sum(y[i_dum:j_dum])
        t1 = x[i_dum]
        t2 = x[j_dum]
        T_div = t2-t1
        return T_div
    
    # index_w = 6
    # index_c = 0
    # # keys_to_plot = [(index_w,0), (index_w,1), (index_w,2), (index_w,3), (index_w,4)]
    # keys_to_plot = [(0,index_c), (1,index_c), (2,index_c), (3,index_c), (4,index_c), (5,index_c), (6,index_c)]
    
    # index_c= 1
    # keys_to_plot = keys_to_plot + [(0,index_c), (1,index_c), (2,index_c), (3,index_c), (4,index_c), (5,index_c), (6,index_c)]
    
    # index_c= 4
    # keys_to_plot = keys_to_plot + [(0,index_c), (1,index_c), (2,index_c), (3,index_c), (4,index_c), (5,index_c), (6,index_c)]
    
    # plt.figure()
    # for key_c in range(len(keys_to_plot)):
    #     key = keys_to_plot[key_c]
    #     lw_ind = key[0]
    #     lc_ind = key[1]
    #     lw = lw_list[lw_ind]
    #     lc = lc_list[lc_ind]
    #     leg = "lw="+str(lw)+"; lc="+str(lc)
        
    #     plt.plot(beta_a_norm_dict[key][0,:], beta_a_norm_dict[key][1,:] , label=leg)
    #     lower = beta_a_norm_dict[key][1,:]-beta_a_norm_dict[key][2,:]
    #     upper = beta_a_norm_dict[key][1,:]+beta_a_norm_dict[key][2,:]
    #     plt.fill_between(beta_a_norm_dict[key][0,:], lower, upper, alpha=0.2)
    
    # plt.plot(beta_a_norm_dict[key][0,:], 1+0.0*beta_a_norm_dict[key][0,:], linestyle='--', color='k')
    # title = "type: "+celltype
    # plt.title(title)
    # plt.xlabel("t (h)")
    # plt.ylabel(r"$\beta_a/\beta_0$")
    # plt.grid()
    # # plt.legend()
    # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # plt.tight_layout()
    # plt.savefig("beta_a_norm.PNG", dpi=300)
    
    if celltype!="C":
        print("Error! use the right function for each type!")
        sys.exit()
    
    labels_list = []
    x_avg_list = []
    y_avg_list = []
    y_err_list = []
    
    for lc_c in range(len(lc_list)):
        lc = lc_list[lc_c]
        labels_list.append("lc="+str(lc))
        
        x_avg_dum = beta_a_norm_dict[(0,lc_c)][0,:]
        
        n_averaging = len(lw_list)
        y_mat_dum = 0.0*np.zeros((n_averaging, len(x_avg_dum)))
        for lw_c in range(len(lw_list)):
            key = (lw_c, lc_c)
            y_mat_dum[lw_c,:] = beta_a_norm_dict[key][1,:]
        y_avg_dum = np.mean(y_mat_dum, axis=0)
        # y_err_dum = np.std(y_mat_dum, axis=0)/np.sqrt(n_averaging-1)
        y_err_dum = np.std(y_mat_dum, axis=0)
        
        x_avg_list.append(x_avg_dum)
        y_avg_list.append(y_avg_dum)
        y_err_list.append(y_err_dum)
    
    for lc_c in range(len(lc_list)):
        lc = lc_list[lc_c]
        data_to_save = np.zeros((3,len(x_avg_list[lc_c])))
        
        data_to_save[0,:] = x_avg_list[lc_c].copy()
        data_to_save[1,:] = y_avg_list[lc_c].copy()
        data_to_save[2,:] = y_err_list[lc_c].copy()
        
        file_name = "beta_norm_"+celltype+"_"+str(lc_c)+".csv"
        np.savetxt("pp_data/"+file_name, X=data_to_save, fmt='%.6e', delimiter=' , ')
    
    plt.figure()
    for lc_c in range(len(lc_list)):
        lc = lc_list[lc_c]
        leg = "lc="+str(lc)
        
        plt.plot(x_avg_list[lc_c], y_avg_list[lc_c] , label=leg)
        lower = y_avg_list[lc_c] - y_err_list[lc_c]
        upper = y_avg_list[lc_c] + y_err_list[lc_c]
        plt.fill_between(x_avg_list[lc_c], lower, upper, alpha=0.2)
    
    plt.plot(x_avg_list[lc_c], 1+0.0*x_avg_list[lc_c], linestyle='--', color='k')
    title = "type: "+celltype
    plt.title(title)
    plt.xlabel("t (h)")
    plt.ylabel(r"$\beta_a/\beta_0$")
    plt.grid()
    # plt.legend()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    plt.savefig("beta_a_norm_C.PNG", dpi=300)
    
    
    # analysis for division time
    T_avg_list = np.zeros(len(lc_list))
    T_err_list = np.zeros(len(lc_list))
    T_vals = np.zeros((2,len(lc_list)))
    beta_0_c = 0.0398
    for lc_c in range(len(lc_list)):
        x_analysis = x_avg_list[lc_c].copy()
        y_avg_analysis = beta_0_c * y_avg_list[lc_c].copy()
        y_upper_analysis = beta_0_c * ( y_avg_list[lc_c].copy()+y_err_list[lc_c].copy() )
        y_lower_analysis = beta_0_c * ( y_avg_list[lc_c].copy()-y_err_list[lc_c].copy() )
        
        dt = x_analysis[1]-x_analysis[0]
        
        T_avg = beta_norm_analysis(x_analysis, y_avg_analysis)
        T_upper = beta_norm_analysis(x_analysis, y_lower_analysis)
        T_lower = beta_norm_analysis(x_analysis, y_upper_analysis)

        T_avg_list[lc_c] = T_avg
        T_err_list[lc_c] = ( ((T_upper-T_lower)/2)**2 + dt**2 ) ** 0.5
    
    T_vals[0,:] = T_avg_list
    T_vals[1,:] = T_err_list
    np.savetxt("pp_data/"+"T_vals.csv", X=T_vals, fmt='%.6e', delimiter=' , ')
    # analysis for division time
    
    
    plt.figure()
    plt.errorbar(lc_list, y=T_avg_list, yerr = T_err_list, color='k'
                 ,markersize=5, fmt='o')
    title = "Max division time in the affected part"
    plt.title(title)
    plt.xscale("log")
    plt.xlabel(r"$\ell_C$ ("+r"$\mu$m)")
    plt.ylabel(r"$T^{\ast}_{a,C}$")
    plt.grid()
    # plt.legend()
    # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    plt.savefig("T_ast_C.PNG", dpi=300)
    
    
    return

def beta_a_norm_plotter_W(beta_a_norm_dict, celltype): # Works only for type W

    
    
    if celltype!="W":
        print("Error! use the right function for each type!")
        sys.exit()
    
    labels_list = []
    x_avg_list = []
    y_avg_list = []
    y_err_list = []
    
    for lw_c in range(len(lw_list)):
        lw = lw_list[lw_c]
        labels_list.append("lw="+str(lw))
        
        x_avg_dum = beta_a_norm_dict[(lw_c,0)][0,:]
        
        n_averaging = len(lc_list)
        y_mat_dum = 0.0*np.zeros((n_averaging, len(x_avg_dum)))
        for lc_c in range(len(lc_list)):
            key = (lw_c, lc_c)
            y_mat_dum[lc_c,:] = beta_a_norm_dict[key][1,:]
        y_avg_dum = np.mean(y_mat_dum, axis=0)
        # y_err_dum = np.std(y_mat_dum, axis=0)/np.sqrt(n_averaging-1)
        y_err_dum = np.std(y_mat_dum, axis=0)
        
        x_avg_list.append(x_avg_dum)
        y_avg_list.append(y_avg_dum)
        y_err_list.append(y_err_dum)
    
    for lw_c in range(len(lw_list)):
        lw = lw_list[lw_c]
        data_to_save = np.zeros((3,len(x_avg_list[lw_c])))
        
        data_to_save[0,:] = x_avg_list[lw_c].copy()
        data_to_save[1,:] = y_avg_list[lw_c].copy()
        data_to_save[2,:] = y_err_list[lw_c].copy()
        
        file_name = "beta_norm_"+celltype+"_"+str(lw_c)+".csv"
        np.savetxt("pp_data/"+file_name, X=data_to_save, fmt='%.6e', delimiter=' , ')
    
    plt.figure()
    for lw_c in range(len(lw_list)):
        lw = lw_list[lw_c]
        leg = "lw="+str(lw)
        
        plt.plot(x_avg_list[lw_c], y_avg_list[lw_c] , label=leg)
        lower = y_avg_list[lw_c] - y_err_list[lw_c]
        upper = y_avg_list[lw_c] + y_err_list[lw_c]
        plt.fill_between(x_avg_list[lw_c], lower, upper, alpha=0.2)
    
    # plt.plot(x_avg_list[lw_c], 1+0.0*x_avg_list[lc_c], linestyle='--', color='k')
    title = "type: "+celltype
    plt.title(title)
    plt.xlabel("t (h)")
    plt.ylabel(r"$\beta_a/\beta_0$")
    plt.grid()
    # plt.legend()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.tight_layout()
    plt.savefig("beta_a_norm_W.PNG", dpi=300)
    
    
    return

def load_object(cls, filename):
    """
    Load attributes from a JSON file and set them on an instance of `cls`.
    Works even if `cls.__init__` takes no arguments.
    """
    with open(filename, "r") as f:
        data = json.load(f)

    obj = cls()  # create instance with default constructor
    for key, value in data.items():
        setattr(obj, key, value)  # assign attributes dynamically
    return obj

class paramsClass:
  # l_w_0 = 0 # how far the effect ON wt cells can reach
  # l_c_0 = 0 # how far the effect ON ca cells can reach
  test = 0

def data_saver(matrix, file_name):
    data_to_save =  np.zeros((len(lw_list)+1, len(lc_list)+1))
    data_to_save[1:,0] = lw_list.copy()
    data_to_save[0,1:] = lc_list.copy()
    data_to_save[1:,1:] = matrix.copy()
    np.savetxt("pp_data/"+file_name, X=data_to_save, fmt='%.5e', delimiter=' , ')
    return


def beta_aff_norm_plotter_both():
    
    plt.figure()
    
    lw_ind_plot_list = [0,1,2,3,4]
    lw_ind_plot_list.append(len(lw_list)-1)
    for lw_ind in lw_ind_plot_list:
        try:
            lw = float(lw_list[lw_ind])
            file_name = "pp_data/"+"beta_norm_W_"+str(lw_ind)+".csv"
            data_load = np.loadtxt(file_name, delimiter=',')
            x_plot = data_load[0,:]
            y_avg  = data_load[1,:]
            y_err  = data_load[2,:]
            y_upper = y_avg+y_err
            y_lower = y_avg-y_err
            label = 'lw='+str(lw)
            plt.plot(x_plot, y_avg, label=label, linestyle='--')
            plt.fill_between(x_plot, y_lower, y_upper, alpha=0.2)
        except:
            continue
    
        
    lc_ind_plot_list = [0,1,2,3]
    lc_ind_plot_list.append(len(lc_list)-1)
    for lc_ind in lc_ind_plot_list:
        try:
            lc = float(lc_list[lc_ind])
            file_name = "pp_data/"+"beta_norm_C_"+str(lc_ind)+".csv"
            data_load = np.loadtxt(file_name, delimiter=',')
            x_plot = data_load[0,:]
            y_avg  = data_load[1,:]
            y_err  = data_load[2,:]
            y_upper = y_avg+y_err
            y_lower = y_avg-y_err
            label = 'lc='+str(lc)
            plt.plot(x_plot, y_avg, label=label)
            plt.fill_between(x_plot, y_lower, y_upper, alpha=0.2)
        except:
            continue
    
    plt.plot(x_plot, 1.0+0*x_plot, linestyle='dashdot', color='k')
    plt.plot(x_plot, 0.0+0*x_plot, linestyle='dashdot', color='k')
    
    # plt.title(title)
    # plt.grid()
    plt.xlabel("t (h)", fontsize=15)
    plt.ylabel(r"$\beta^a(t)/\beta^0$", fontsize=15)
    # plt.grid()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.savefig("beta_aff_norm_both.PNG", dpi=300)
    # plt.close()
    
    return

lw_list = []
lc_list = []

with open("lw_lc_lists.sh", "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith("lw_list"):
            values = line.split("=", 1)[1].strip("() ;")
            lw_list = [float(x) for x in values.split()]
        elif line.startswith("lc_list"):
            values = line.split("=", 1)[1].strip("() ;")
            lc_list = [float(x) for x in values.split()]



#n_samples = 10
n_samples = int(np.loadtxt("n_samples.txt", delimiter=','))


cost_data_avg = np.zeros((len(lw_list), len(lc_list)))
cost_data_err = np.zeros((len(lw_list), len(lc_list)))

main_cost_data_w_avg = np.zeros((len(lw_list), len(lc_list)))
main_cost_data_c_avg = np.zeros((len(lw_list), len(lc_list)))
main_cost_data_w_err = np.zeros((len(lw_list), len(lc_list)))
main_cost_data_c_err = np.zeros((len(lw_list), len(lc_list)))

deriv_cost_data_w_avg = np.zeros((len(lw_list), len(lc_list)))
deriv_cost_data_c_avg = np.zeros((len(lw_list), len(lc_list)))
deriv_cost_data_w_err = np.zeros((len(lw_list), len(lc_list)))
deriv_cost_data_c_err = np.zeros((len(lw_list), len(lc_list)))


b_w_data_avg = np.zeros((len(lw_list), len(lc_list)))
b_w_data_err = np.zeros((len(lw_list), len(lc_list)))

b_c_data_avg = np.zeros((len(lw_list), len(lc_list)))
b_c_data_err = np.zeros((len(lw_list), len(lc_list)))

# tau_w_data_avg = np.zeros((len(lw_list), len(lc_list)))
# tau_w_data_err = np.zeros((len(lw_list), len(lc_list)))

# tau_c_data_avg = np.zeros((len(lw_list), len(lc_list)))
# tau_c_data_err = np.zeros((len(lw_list), len(lc_list)))



## beta_a ratios
beta_w_a_norm_dict = dict()
beta_c_a_norm_dict = dict()
## beta_a ratios

W_v_to_C_a_dict = dict()


# exp data
overal_WT_pure = np.loadtxt("exp_data/"+"overal_WT_pure.csv", delimiter=',')
overal_WT_mix = np.loadtxt("exp_data/"+"overal_WT_mix.csv", delimiter=',')
overal_C_pure = np.loadtxt("exp_data/"+"overal_C_pure.csv", delimiter=',')
overal_C_mix = np.loadtxt("exp_data/"+"overal_C_mix.csv", delimiter=',')

time_data   = overal_WT_pure[0,:]
WT_norm_avg = overal_WT_pure[1,:]
WT_norm_err = overal_WT_pure[2,:]
C_norm_avg = overal_C_pure[1,:]
C_norm_err = overal_C_pure[2,:]

WT_mix_norm_avg = overal_WT_mix[1,:]
WT_mix_norm_err = overal_WT_mix[2,:]
C_mix_norm_avg  = overal_C_mix[1,:]
C_mix_norm_err  = overal_C_mix[2,:]
# exp data

W_mix_fit_smples_coefs = np.loadtxt("W_mix_fit_smples_coefs.csv", delimiter=',')
C_mix_fit_smples_coefs = np.loadtxt("C_mix_fit_smples_coefs.csv", delimiter=',')

for lw_ind in range(len(lw_list)):
    lw = float(lw_list[lw_ind])
    
    for lc_ind in range(len(lc_list)):
        lc = float(lc_list[lc_ind])
        
        
        folder_name = "lw_"+str(lw)+"__lc_"+str(lc)
        
        params = load_object(paramsClass, folder_name+"/sample_1/params.txt")

        
        time = np.loadtxt(folder_name+"/sample_1/data/time.txt", delimiter=',')
        dt = time[1]-time[0]
        
        A_w_stack = np.empty((0,len(time)))
        A_c_stack = np.empty((0,len(time)))
        
        A_w_v_stack = np.empty((0,len(time)))
        A_c_a_stack = np.empty((0,len(time)))
        
        
        cost_list = [] # total cost
        cost_list.clear()
        cost_list = []
        
        # cost details
        main_cost_list_w = []
        main_cost_list_w.clear()
        main_cost_list_w = []
        
        deriv_cost_list_w = []
        deriv_cost_list_w.clear()
        deriv_cost_list_w = []
        
        main_cost_list_c = []
        main_cost_list_c.clear()
        main_cost_list_c = []
        
        deriv_cost_list_c = []
        deriv_cost_list_c.clear()
        deriv_cost_list_c = []
        # cost details
        
        
        
        b_w_list = []
        b_w_list.clear()
        b_w_list = []
        
        b_c_list = []
        b_c_list.clear()
        b_c_list = []
        
        # tau_w_list = []
        # tau_w_list.clear()
        # tau_w_list = []
        
        # tau_c_list = []
        # tau_c_list.clear()
        # tau_c_list = []
        
        dict_key = tuple((lw_ind, lc_ind))
        beta_w_a_norm_dict[dict_key] = np.zeros((3, len(time)))
        beta_c_a_norm_dict[dict_key] = np.zeros((3, len(time)))
        beta_w_a_norm_dict[dict_key][0,:] = time.copy()
        beta_c_a_norm_dict[dict_key][0,:] = time.copy()
        
        W_v_to_C_a_dict[dict_key] = np.zeros((3, len(time)))
        W_v_to_C_a_dict[dict_key][0,:] = time.copy()
        
        beta_w_a_norm_dummy_all_samples = 0.0*np.zeros((n_samples, len(time)))
        beta_c_a_norm_dummy_all_samples = 0.0*np.zeros((n_samples, len(time)))
        
        for sample_c in range(n_samples):
            address = folder_name + "/sample_"+str(sample_c+1)
            GD_log_data = np.loadtxt(address+"/GD_log.csv", delimiter=',')
            cost = GD_log_data[-1, 0]
            cost_list.append(cost)
            
            #reading details of cost
            with open(address+"/cost_log.txt", 'r') as f:
                content = f.read()
            
            # Split into blocks using the separator
            blocks = content.strip().split('--------')
            
            # Take the last non-empty block
            last_block = blocks[-1].strip()
            if last_block == '':
                last_block = blocks[-2].strip()
            
            # Split lines and convert to float
            cost_detail_data = []
            for line in last_block.splitlines():
                if line.strip():  # skip empty lines
                    numbers = list(map(float, line.split()))
                    cost_detail_data.append(numbers)
            cost_detail_data=np.array(cost_detail_data)
            
            main_cost_list_w.append(cost_detail_data[0,0])
            main_cost_list_c.append(cost_detail_data[1,0])
            deriv_cost_list_w.append(cost_detail_data[0,1])
            deriv_cost_list_c.append(cost_detail_data[1,1])
            #reading details of cost
            
            A_w_sample = np.loadtxt(address+"/data/"+"A_w_mat.txt", delimiter=',')
            A_w_stack = np.vstack([A_w_stack, A_w_sample])
            
            A_c_sample = np.loadtxt(address+"/data/"+"A_c_mat.txt", delimiter=',')
            A_c_stack = np.vstack([A_c_stack, A_c_sample])
            
            A_w_v_sample = np.loadtxt(address+"/data/"+"A_w_v_mat.txt", delimiter=',')
            A_w_v_stack = np.vstack([A_w_v_stack, A_w_v_sample])
            
            A_c_a_sample = np.loadtxt(address+"/data/"+"A_c_aff_mat.txt", delimiter=',')
            A_c_a_stack = np.vstack([A_c_a_stack, A_c_a_sample])
            
            
            # print(folder_name)
            
            beta_w_aff_mat = np.loadtxt(address+"/data/"+'beta_w_aff_mat.txt',  delimiter=',')
            beta_c_aff_mat = np.loadtxt(address+"/data/"+'beta_c_aff_mat.txt',  delimiter=',')
            
            beta_w_aff_norm = beta_w_aff_mat / 0.0284
            beta_c_aff_norm = beta_c_aff_mat / 0.0398
            
            beta_w_aff_norm_avg_sample = np.mean(beta_w_aff_norm, axis=0)
            beta_c_aff_norm_avg_sample = np.mean(beta_c_aff_norm, axis=0)
            
            beta_w_a_norm_dummy_all_samples[sample_c, :] = beta_w_aff_norm_avg_sample.copy()
            beta_c_a_norm_dummy_all_samples[sample_c, :] = beta_c_aff_norm_avg_sample.copy()
            
            # # cost for first deriv difference
            # A_w_norm_sample = A_w_sample / A_w_sample[:, [0]]
            # A_w_norm_sample_avg = np.mean(A_w_norm_sample, axis=0)
            # spline_w = UnivariateSpline(time, np.log(A_w_norm_sample_avg), s=0)
            # f1_w = spline_w.derivative(1)(time)
            # f2_w = spline_w.derivative(2)(time)
            # curvature_w = f2_w / (1 + f1_w**2)**1.5
            
            # fit_coef_sample = np.random.randint(np.shape(W_mix_fit_smples_coefs)[1])
            # a = W_mix_fit_smples_coefs[0,fit_coef_sample]
            # b = W_mix_fit_smples_coefs[1,fit_coef_sample]
            # c = W_mix_fit_smples_coefs[2,fit_coef_sample]
            
            # f1_w_exp_fit = 2*a*time+b
            # deriv_cost = sum(dt*(f1_w-f1_w_exp_fit)**2)
            # deriv_cost_list.append(deriv_cost)
            # # cost for first deriv difference
            
            # # cost for curv
            # fit_coef_sample = np.random.randint(np.shape(W_mix_fit_smples_coefs)[1])
            # a = W_mix_fit_smples_coefs[0,fit_coef_sample]
            # b = W_mix_fit_smples_coefs[1,fit_coef_sample]
            # c = W_mix_fit_smples_coefs[2,fit_coef_sample]
            
            # f2_w_exp_fit = 2*a
            # fit_curv_w = f2_w_exp_fit / (1+f1_w_exp_fit**2)**1.5
            # curv_cost = sum(dt*(curvature_w-fit_curv_w)**2)
            # curv_cost_list.append(curv_cost)
            # # cost for curv
            
            # b_w, b_c, tau_w, tau_c
            b_w = GD_log_data[-1,4]
            b_c = GD_log_data[-1,5]
            # tau_w = GD_log_data[-1,5]
            # tau_c = GD_log_data[-1,6]
            b_w_list.append(b_w)
            b_c_list.append(b_c)
            # tau_w_list.append(tau_w)
            # tau_c_list.append(tau_c)
            # b_w, b_c, tau_w, tau_c
            
        
        # # plot b_w, b_c histograms
        # plt.figure()
        # plt.hist(b_w_list)
        # plt.title('b_w')
        
        # plt.figure()
        # plt.hist(b_c_list)
        # plt.title('b_c')
        # # plot b_w, b_c histograms
        
        A_w_norm_stack = A_w_stack / A_w_stack[:, [0]]
        A_c_norm_stack = A_c_stack / A_c_stack[:, [0]]
        
        A_w_norm_avg = np.mean(A_w_norm_stack, axis=0)
        A_w_norm_err = np.std(A_w_norm_stack, axis=0)/np.sqrt(np.shape(A_w_norm_stack)[1])
        
        A_c_norm_avg = np.mean(A_c_norm_stack, axis=0)
        A_c_norm_err = np.std(A_c_norm_stack, axis=0)/np.sqrt(np.shape(A_c_norm_stack)[1])
        
        beta_w_a_norm_dict[dict_key][1,:] = np.mean(beta_w_a_norm_dummy_all_samples, axis=0)
        beta_w_a_norm_dict[dict_key][2,:] = np.std(beta_w_a_norm_dummy_all_samples, axis=0)/np.sqrt(n_samples-1)
        beta_c_a_norm_dict[dict_key][1,:] = np.mean(beta_c_a_norm_dummy_all_samples, axis=0)
        beta_c_a_norm_dict[dict_key][2,:] = np.std(beta_c_a_norm_dummy_all_samples, axis=0)/np.sqrt(n_samples-1)
        
        
        W_v_to_C_a_stack = (A_w_v_stack/params.a_w)/(A_c_a_stack/params.a_c)
        W_v_to_C_a_dict[dict_key][1,:] = np.mean(W_v_to_C_a_stack, axis=0)
        W_v_to_C_a_dict[dict_key][2,:] = np.std(W_v_to_C_a_stack, axis=0)
        
        plt.figure()
        title = folder_name
        x_plot = W_v_to_C_a_dict[dict_key][0,:]
        y_avg = np.mean(W_v_to_C_a_stack, axis=0)
        y_err = np.std(W_v_to_C_a_stack, axis=0)
        y_upper = y_avg+y_err
        y_lower = y_avg-y_err
        plt.plot(x_plot, y_avg)
        plt.fill_between(x_plot, y_lower, y_upper, alpha=0.2)
        plt.title(title)
        plt.xlabel("t (h)")
        plt.ylabel(r"$W^v(t)/C^a(t)$")
        plt.grid()
        # plt.legend()
        plt.tight_layout()
        plt.savefig(folder_name+"/"+title+".PNG", dpi=300)
        plt.close()
        
        
        # ## plot for this lw , lc
        # plt.figure()
        # err1 = plt.errorbar(time_data, WT_norm_avg,  yerr=WT_norm_err, fmt='o', color='m', ecolor='m', capsize=2, label='pure WT')
        # err2 =plt.errorbar(time_data, C_norm_avg,  yerr=C_norm_err, fmt='o', color='g', ecolor='g', capsize=2, label='pure C')
        # err3 = plt.errorbar(time_data, WT_mix_norm_avg, yerr=WT_mix_norm_err, fmt='s', color='m', ecolor='m', capsize=2, label='mixed WT', mfc='none')
        # err4 =plt.errorbar(time_data, C_mix_norm_avg, yerr=C_mix_norm_err, fmt='s', color='g', ecolor='g', capsize=2, label='mixed C', mfc='none')
        # plt.errorbar(time, y=A_w_norm_avg, yerr = A_w_norm_err, label='w', color='m')
        # plt.errorbar(time, y=A_c_norm_avg, yerr = A_c_norm_err, label='c', color='g')
        # plt.title(folder_name)
        # plt.yscale("log")
        # plt.grid()
        # plt.legend(fontsize=15)
        # plt.xlabel("time", fontsize=15)
        # plt.ylabel("N(t)/N(0)", fontsize=15)
        # ## plot for this lw , lc
        
        
        # # curvature
        # spline = UnivariateSpline(time, np.log(A_w_norm_avg), s=0)
        # # First and second derivatives
        # f1 = spline.derivative(1)(time)
        # f2 = spline.derivative(2)(time)
        # curvature = f2 / (1 + f1**2)**1.5
        
        # n_discrete = 1000
        # x_smooth = np.linspace(time.min(), time.max(), n_discrete)
        # y_smooth = spline(x_smooth)

        # W_mix_fit_smples_coefs = np.loadtxt("W_mix_fit_smples_coefs.csv", delimiter=',')
        # fit_coef_samples = np.random.randint(0,np.shape(W_mix_fit_smples_coefs)[1],n_samples)
        
        # plt.figure()
        # for i in range(n_samples):
        #     a = W_mix_fit_smples_coefs[0,fit_coef_samples[i]]
        #     b = W_mix_fit_smples_coefs[1,fit_coef_samples[i]]
        #     c = W_mix_fit_smples_coefs[2,fit_coef_samples[i]]
            
        #     fit_prime = 2 * a* x_smooth + b
        #     fit_d_prime = 2*a
        #     fit_curv = fit_d_prime / (1+fit_prime**2)**1.5
            
        #     plt.plot(x_smooth, a*x_smooth**2+b*x_smooth+c, linestyle='--', color='r')
            
        # plt.errorbar(time, y=np.log(A_w_norm_avg), yerr = A_w_norm_err/A_w_norm_avg, label='w', color='m')
        # plt.plot(x_smooth, y_smooth, '-', label='Spline fit', linewidth=2, zorder=10, alpha=0.5)  # spline curve
        # plt.xlabel("t")
        # # plt.ylabel("W(t)/W(0)")
        # plt.grid(True)
        # plt.legend()
        # plt.show()
        
        
        # plt.figure()
        # for i in range(n_samples):
        #     a = W_mix_fit_smples_coefs[0,fit_coef_samples[i]]
        #     b = W_mix_fit_smples_coefs[1,fit_coef_samples[i]]
        #     c = W_mix_fit_smples_coefs[2,fit_coef_samples[i]]
            
        #     fit_prime = 2 * a* x_smooth + b
            
        #     plt.plot(x_smooth, 2*a*x_smooth+b, linestyle='--', color='r')
            
        # plt.plot(time, f1, '-', linewidth=2, zorder=10, alpha=0.5)  # spline curve
        # plt.xlabel("t")
        # plt.ylabel("W(t)/W(0)")
        # plt.grid(True)
        # plt.legend()
        # plt.show()
        
        
        # plt.figure()
        # for i in range(n_samples):
        #     a = W_mix_fit_smples_coefs[0,fit_coef_samples[i]]
        #     b = W_mix_fit_smples_coefs[1,fit_coef_samples[i]]
        #     c = W_mix_fit_smples_coefs[2,fit_coef_samples[i]]
            
        #     fit_prime = 2 * a* x_smooth + b
        #     fit_d_prime = 2*a
        #     fit_curv = fit_d_prime / (1+fit_prime**2)**1.5
            
        #     plt.plot(x_smooth, fit_curv, linestyle='--', color='r')
            
        # plt.plot(time, curvature, '-', label='Spline fit', linewidth=2, zorder=10, alpha=0.5)  # spline curve
        # plt.xlabel("t")
        # plt.ylabel("curv_W")
        # plt.grid(True)
        # plt.legend()
        # plt.show()
        # # curvature
        
        cost_avg = np.mean(cost_list)
        cost_err = np.std(cost_list)/np.sqrt(n_samples-1)
        cost_data_avg[lw_ind, lc_ind] = cost_avg
        cost_data_err[lw_ind, lc_ind] = cost_err
        
        main_cost_data_w_avg[lw_ind, lc_ind] = np.mean(main_cost_list_w)
        main_cost_data_w_err[lw_ind, lc_ind] = np.std(main_cost_list_w)/np.sqrt(n_samples-1)
        deriv_cost_data_w_avg[lw_ind, lc_ind] = np.mean(deriv_cost_list_w)
        deriv_cost_data_w_err[lw_ind, lc_ind] = np.std(deriv_cost_list_w)/np.sqrt(n_samples-1)
        
        main_cost_data_c_avg[lw_ind, lc_ind] = np.mean(main_cost_list_c)
        main_cost_data_c_err[lw_ind, lc_ind] = np.std(main_cost_list_c)/np.sqrt(n_samples-1)
        deriv_cost_data_c_avg[lw_ind, lc_ind] = np.mean(deriv_cost_list_c)
        deriv_cost_data_c_err[lw_ind, lc_ind] = np.std(deriv_cost_list_c)/np.sqrt(n_samples-1)
        
        # deriv_cost_avg = np.mean(deriv_cost_list)
        # deriv_cost_err = np.std(deriv_cost_list)/np.sqrt(n_samples-1)
        # deriv_cost_data_avg[lw_ind, lc_ind] = deriv_cost_avg
        # deriv_cost_data_err[lw_ind, lc_ind] = deriv_cost_err
        
        # curv_cost_avg = np.mean(curv_cost_list)
        # curv_cost_err = np.std(curv_cost_list)/np.sqrt(n_samples-1)
        # curv_cost_data_avg[lw_ind, lc_ind] = curv_cost_avg
        # curv_cost_data_err[lw_ind, lc_ind] = curv_cost_err
        
        b_w_avg = np.mean(b_w_list)
        b_w_err = np.std(b_w_list)/np.sqrt(n_samples-1)
        b_w_data_avg[lw_ind, lc_ind] = b_w_avg
        b_w_data_err[lw_ind, lc_ind] = b_w_err
        
        b_c_avg = np.mean(b_c_list)
        b_c_err = np.std(b_c_list)/np.sqrt(n_samples-1)
        b_c_data_avg[lw_ind, lc_ind] = b_c_avg
        b_c_data_err[lw_ind, lc_ind] = b_c_err
        
        # tau_w_avg = np.mean(tau_w_list)
        # tau_w_err = np.std(tau_w_list)/np.sqrt(n_samples-1)
        # tau_w_data_avg[lw_ind, lc_ind] = tau_w_avg
        # tau_w_data_err[lw_ind, lc_ind] = tau_w_err
        
        # tau_c_avg = np.mean(tau_c_list)
        # tau_c_err = np.std(tau_c_list)/np.sqrt(n_samples-1)
        # tau_c_data_avg[lw_ind, lc_ind] = tau_c_avg
        # tau_c_data_err[lw_ind, lc_ind] = tau_c_err


# plot W_v(t)/C_a(t) for sample lc, lw
try:
    selected_keys = [(3,0), (3,1), (3,4)] #lw=30; lc=11,33,1100
    plt.figure()
    title = 'all for lw=30'
    for key in selected_keys:
        
        x_plot = W_v_to_C_a_dict[key][0,:]
        y_avg = W_v_to_C_a_dict[key][1,:]
        y_err = W_v_to_C_a_dict[key][2,:]
        y_upper = y_avg+y_err
        y_lower = y_avg-y_err
        label = 'lc='+str(lc_list[key[1]])
        plt.plot(x_plot, y_avg, label = label)
        plt.fill_between(x_plot, y_lower, y_upper, alpha=0.2)
    plt.title(title)
    plt.xlabel("t (h)")
    plt.ylabel(r"$W^v(t)/C^a(t)$")
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.savefig("W_v_to_C_a.PNG", dpi=300)
    plt.close()
except:
    pass
# plot W_v(t)/C_a(t) for sample lc, lw

b_w_difference_avg = np.zeros((len(lw_list), len(lc_list)))
b_w_difference_err = np.zeros((len(lw_list), len(lc_list)))

b_c_difference_avg = np.zeros((len(lw_list), len(lc_list)))
b_c_difference_err = np.zeros((len(lw_list), len(lc_list)))

b_w_ratio_avg = np.zeros((len(lw_list), len(lc_list)))
b_w_ratio_err = np.zeros((len(lw_list), len(lc_list)))

b_c_ratio_avg = np.zeros((len(lw_list), len(lc_list)))
b_c_ratio_err = np.zeros((len(lw_list), len(lc_list)))



for lw_ind in range(len(lw_list)):
    for lc_ind in range(len(lc_list)):
        
        b_w_difference_avg[lw_ind, lc_ind] = b_w_data_avg[lw_ind, lc_ind]-b_w_data_avg[-1, lc_ind]
        b_w_difference_err[lw_ind, lc_ind] = b_w_data_err[lw_ind, lc_ind]+b_w_data_err[-1, lc_ind]

        b_c_difference_avg[lw_ind, lc_ind] = b_c_data_avg[lw_ind, lc_ind]-b_c_data_avg[lw_ind, -1]
        b_c_difference_err[lw_ind, lc_ind] = b_c_data_err[lw_ind, lc_ind]+b_c_data_err[lw_ind, -1]
        

        b_w_ratio_avg[lw_ind, lc_ind] = b_w_data_avg[lw_ind, lc_ind]/b_w_data_avg[-1, lc_ind]
        b_w_ratio_err[lw_ind, lc_ind] = (b_w_data_avg[lw_ind, lc_ind]*b_w_data_err[-1, lc_ind] + \
                                         b_w_data_err[lw_ind, lc_ind]*b_w_data_avg[-1, lc_ind]) / \
                                         b_w_data_avg[-1, lc_ind] **2 ;
            

        b_c_ratio_avg[lw_ind, lc_ind] = b_c_data_avg[lw_ind, lc_ind]/b_c_data_avg[lw_ind, -1]
        b_c_ratio_err[lw_ind, lc_ind] = (b_c_data_avg[lw_ind, lc_ind]*b_c_data_err[lw_ind, -1] + \
                                         b_c_data_err[lw_ind, lc_ind]*b_c_data_avg[lw_ind, -1]) / \
                                         b_c_data_avg[lw_ind, -1] **2 ;
        
        
        
        

try:
    os.makedirs("pp_data")
except:
    pass
    
data_saver(cost_data_avg, "cost_data_avg.csv")
data_saver(cost_data_err, "cost_data_err.csv")

data_saver(main_cost_data_w_avg, "main_cost_data_w_avg.csv")
data_saver(main_cost_data_w_err, "main_cost_data_w_err.csv")
data_saver(deriv_cost_data_w_avg, "deriv_cost_data_w_avg.csv")
data_saver(deriv_cost_data_w_err, "deriv_cost_data_w_err.csv")

data_saver(main_cost_data_c_avg, "main_cost_data_c_avg.csv")
data_saver(main_cost_data_c_err, "main_cost_data_c_err.csv")
data_saver(deriv_cost_data_c_avg, "deriv_cost_data_c_avg.csv")
data_saver(deriv_cost_data_c_err, "deriv_cost_data_c_err.csv")

data_saver(b_w_data_avg, "b_w_data_avg.csv")
data_saver(b_w_data_err, "b_w_data_err.csv")

data_saver(b_c_data_avg, "b_c_data_avg.csv")
data_saver(b_c_data_err, "b_c_data_err.csv")

# data_saver(tau_w_data_avg, "tau_w_data_avg.csv")
# data_saver(tau_w_data_err, "tau_w_data_err.csv")

# data_saver(tau_c_data_avg, "tau_c_data_avg.csv")
# data_saver(tau_c_data_err, "tau_c_data_err.csv")

np.savetxt("pp_data/"+"lw_list.csv", X=lw_list, fmt='%.5e', delimiter=' , ')
np.savetxt("pp_data/"+"lc_list.csv", X=lc_list, fmt='%.5e', delimiter=' , ')
b_w_tot_avg = np.zeros(len(lw_list))
b_w_tot_err = np.zeros(len(lw_list))
b_c_tot_avg = np.zeros(len(lc_list))
b_c_tot_err = np.zeros(len(lc_list))

for lw_ind in range(len(lw_list)):
    b_w_tot_avg[lw_ind] = np.mean(b_w_data_avg[lw_ind,:])
    b_w_tot_err[lw_ind] = np.std(b_w_data_avg[lw_ind,:])/np.sqrt(len(lc_list)-1)
    
for lc_ind in range(len(lc_list)):
    b_c_tot_avg[lc_ind] = np.mean(b_c_data_avg[:,lc_ind])
    b_c_tot_err[lc_ind] = np.std(b_c_data_avg[:,lc_ind])/np.sqrt(len(lw_list)-1)

np.savetxt("pp_data/"+"b_w_tot_avg.csv", X=b_w_tot_avg, fmt='%.5e', delimiter=' , ')
np.savetxt("pp_data/"+"b_w_tot_err.csv", X=b_w_tot_err, fmt='%.5e', delimiter=' , ')
np.savetxt("pp_data/"+"b_c_tot_avg.csv", X=b_c_tot_avg, fmt='%.5e', delimiter=' , ')
np.savetxt("pp_data/"+"b_c_tot_err.csv", X=b_c_tot_err, fmt='%.5e', delimiter=' , ')



beta_a_norm_plotter_C(beta_c_a_norm_dict, "C")
beta_a_norm_plotter_W(beta_w_a_norm_dict, "W")

# sdiuysgd

# # plot b_w_tot
# plt.figure()
# # plt.errorbar(lw_list, y=b_w_tot_avg-b_w_tot_avg[-1], yerr = b_w_tot_err, label='tot_avg', color='k',markersize=5, fmt='o')
# plt.errorbar(lw_list, y=b_w_tot_avg, yerr = b_w_tot_err, label='tot_avg', color='k',markersize=5, fmt='o')
# plt.xlabel("l_w", fontsize=15)
# plt.ylabel(r"$\overline{b_w}$", fontsize=15)
# # plt.yscale("log")
# plt.xscale("log")
# plt.yscale("log")
# plt.grid()
# plt.tight_layout()
# plt.savefig("b_w_tot_avg.PNG", dpi=300)
# # plot b_w_tot

# # plot b_c_tot
# plt.figure()
# # plt.errorbar(lc_list, y=b_c_tot_avg-b_c_tot_avg[-1], yerr = b_c_tot_err, label='tot_avg', color='k',markersize=5, fmt='o')
# plt.errorbar(lc_list, y=b_c_tot_avg, yerr = b_c_tot_err, label='tot_avg', color='k',markersize=5, fmt='o')
# plt.xlabel("l_c", fontsize=15)
# plt.ylabel(r"$\overline{b_c}$", fontsize=15)
# # plt.yscale("log")
# plt.xscale("log")
# plt.yscale("log")
# plt.grid()
# plt.tight_layout()
# plt.savefig("b_c_tot_avg.PNG", dpi=300)
# # plot b_c_tot


# plot cost
plt.figure()
for i in range(len(lc_list)):
    label = 'lc = '+str(lc_list[i])
    plt.errorbar(lw_list, y=cost_data_avg[:,i], yerr = cost_data_err[:,i], label=label)
plt.legend(fontsize=15)
plt.xlabel("l_w", fontsize=15)
plt.ylabel("cost", fontsize=15)
plt.xscale("log")
plt.grid()
plt.tight_layout()
plt.savefig("cost.PNG", dpi=300)
# plot cost

# plot cost components
fig, axs = plt.subplots(2, 2, figsize=(8, 6))
# Top-left
for i in range(len(lc_list)):
    label = 'lc = '+str(lc_list[i])
    axs[0, 0].errorbar(lw_list, y=main_cost_data_w_avg[:,i], yerr = main_cost_data_w_err[:,i], label=label)
axs[0, 0].set_title("main, w")
axs[0, 0].set_xlabel("l_w", fontsize=15)
axs[0, 0].set_ylabel("cost", fontsize=15)
axs[0, 0].set_xscale("log")
axs[0, 0].grid()
# Top-right
for i in range(len(lc_list)):
    label = 'lc = '+str(lc_list[i])
    axs[0, 1].errorbar(lw_list, y=main_cost_data_c_avg[:,i], yerr = main_cost_data_c_err[:,i], label=label)
axs[0, 1].set_title("main, c")
axs[0, 1].set_xlabel("l_w", fontsize=15)
axs[0, 1].set_ylabel("cost", fontsize=15)
axs[0, 1].set_xscale("log")
axs[0, 1].grid()
# Bottom-left
for i in range(len(lc_list)):
    label = 'lc = '+str(lc_list[i])
    axs[1, 0].errorbar(lw_list, y=deriv_cost_data_w_avg[:,i], yerr = deriv_cost_data_w_err[:,i], label=label)
axs[1, 0].set_title("deriv, w")
axs[1, 0].set_xlabel("l_w", fontsize=15)
axs[1, 0].set_ylabel("cost", fontsize=15)
axs[1, 0].set_xscale("log")
axs[1, 0].grid()
# Bottom-right
for i in range(len(lc_list)):
    label = 'lc = '+str(lc_list[i])
    axs[1, 1].errorbar(lw_list, y=deriv_cost_data_c_avg[:,i], yerr = deriv_cost_data_c_err[:,i], label=label)
axs[1, 1].set_title("deriv, c")
axs[1, 1].set_xlabel("l_w", fontsize=15)
axs[1, 1].set_ylabel("cost", fontsize=15)
axs[1, 1].set_xscale("log")
axs[1, 1].grid()
# Adjust layout
plt.tight_layout()
plt.legend()
plt.savefig("cost_details.PNG", dpi=400)
# plot cost components


# plot cost avg
plt.figure()
max_vals = np.max(cost_data_avg + cost_data_err, axis=1)
min_vals = np.min(cost_data_avg - cost_data_err, axis=1)
err_vals = (max_vals - min_vals) / 2
# err_vals = np.mean(b_w_difference_err, axis=1)
plt.errorbar(lw_list, y=np.mean(cost_data_avg, axis=1), yerr = err_vals, label='average', color='k',markersize=5, fmt='o')
plt.plot(lw_list, np.mean(cost_data_avg, axis=1), color='k', linestyle='--')
# log_x = np.log10(lw_list)
# log_y = np.log10(np.mean(b_w_difference_avg, axis=1))
# log_sig = err_vals / (np.log(10) * np.mean(b_w_difference_avg, axis=1)) 
# w = 1 / log_sig**2
# max_index = 7
# coef = np.polyfit(log_x[0:max_index], log_y[0:max_index], deg=1, w=w[0:max_index])   # highest power first
# lw_fit_plot = np.linspace(lw_list[0], lw_list[max_index], 1000)
# bw_diff_fit_plot = (lw_fit_plot**coef[0])*(10**coef[1])
# plt.plot(lw_fit_plot, bw_diff_fit_plot, color='r', linestyle='--')
plt.legend(fontsize=10)
plt.xlabel("l_w", fontsize=15)
plt.ylabel("cost (avg)", fontsize=15)
# plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.tight_layout()
plt.savefig("cost_avg.PNG", dpi=300)
# plot cost avg


# # plot cost
# plt.figure()
# for i in range(len(lc_list)):
#     label = 'lc = '+str(lc_list[i])
#     plt.errorbar(lw_list, y=deriv_cost_data_avg[:,i], yerr = deriv_cost_data_err[:,i], label=label)
# plt.legend(fontsize=15)
# plt.xlabel("l_w", fontsize=15)
# plt.ylabel("deriv_cost_w", fontsize=15)
# plt.xscale("log")
# plt.grid()
# plt.tight_layout()
# plt.savefig("deriv_cost_w.PNG", dpi=300)
# # plot cost


# # plot cost
# plt.figure()
# for i in range(len(lc_list)):
#     label = 'lc = '+str(lc_list[i])
#     plt.errorbar(lw_list, y=curv_cost_data_avg[:,i], yerr = curv_cost_data_err[:,i], label=label)
# plt.legend(fontsize=15)
# plt.xlabel("l_w", fontsize=15)
# plt.ylabel("curv_cost_w", fontsize=15)
# plt.xscale("log")
# plt.grid()
# plt.tight_layout()
# plt.savefig("curv_cost_w.PNG", dpi=300)
# # plot cost


# plot b_w
plt.figure()
for i in range(len(lc_list)):
    label = 'lc = '+str(lc_list[i])
    plt.errorbar(lw_list, y=b_w_data_avg[:,i], yerr = b_w_data_err[:,i], label=label)
plt.legend(fontsize=10)
plt.xlabel("l_w", fontsize=15)
plt.ylabel("b_w", fontsize=15)
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.tight_layout()
plt.savefig("b_w.PNG", dpi=300)
# plot b_w

# plot b_w-b_w_terminal
plt.figure()
for i in range(len(lc_list)):
    label = 'lc = '+str(lc_list[i])
    plt.errorbar(lw_list, y=b_w_difference_avg[:,i], yerr = b_w_difference_err[:,i], label=label)
plt.legend(fontsize=10)
plt.xlabel("l_w", fontsize=15)
plt.ylabel("b_w_difference", fontsize=15)
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.tight_layout()
plt.savefig("b_w_difference.PNG", dpi=300)
# plot b_w-b_w_terminal




# plot b_w-b_w_terminal (avg)
plt.figure()
max_vals = np.max(b_w_difference_avg + b_w_difference_err, axis=1)
min_vals = np.min(b_w_difference_avg - b_w_difference_err, axis=1)
err_vals = (max_vals - min_vals) / 2
# err_vals = np.mean(b_w_difference_err, axis=1)
plt.errorbar(lw_list, y=np.mean(b_w_difference_avg, axis=1), yerr = err_vals, label='average', color='k',markersize=5, fmt='o')
log_x = np.log10(lw_list)
log_y = np.log10(np.mean(b_w_difference_avg, axis=1))
log_sig = err_vals / (np.log(10) * np.mean(b_w_difference_avg, axis=1)) 
w = 1 / log_sig**2
max_index = 7
try:
    coef = np.polyfit(log_x[0:max_index], log_y[0:max_index], deg=1, w=w[0:max_index])   # highest power first
    lw_fit_plot = np.linspace(lw_list[0], lw_list[max_index], 1000)
    bw_diff_fit_plot = (lw_fit_plot**coef[0])*(10**coef[1])
    plt.plot(lw_fit_plot, bw_diff_fit_plot, color='r', linestyle='--')
except:
    pass
plt.legend(fontsize=10)
plt.xlabel("l_w", fontsize=15)
plt.ylabel("b_w_difference (avg)", fontsize=15)
plt.yscale("log")
plt.xscale("log")
plt.title("slope: "+str(round(coef[0],2)))
plt.grid()
plt.tight_layout()
plt.savefig("b_w_difference_avg.PNG", dpi=300)
# plot b_w-b_w_terminal (avg)


# plot b_w-b_w_terminal
plt.figure()
for i in range(len(lc_list)):
    label = 'lc = '+str(lc_list[i])
    plt.errorbar(lw_list, y=b_w_ratio_avg[:,i], yerr = b_w_ratio_err[:,i], label=label)
plt.legend(fontsize=10)
plt.xlabel("l_w", fontsize=15)
plt.ylabel("b_w_ratio", fontsize=15)
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.tight_layout()
plt.savefig("b_w_ratio.PNG", dpi=300)
# plot b_w-b_w_terminal



# # plot b_c
# plt.figure()
# for i in range(len(lc_list)):
#     label = 'lc = '+str(lc_list[i])
#     plt.errorbar(lw_list, y=b_c_data_avg[:,i], yerr = b_c_data_err[:,i], label=label)
# plt.legend(fontsize=15)
# plt.xlabel("l_w", fontsize=15)
# plt.ylabel("b_c", fontsize=15)
# plt.xscale("log")
# plt.grid()
# plt.tight_layout()
# plt.savefig("b_c.PNG", dpi=300)
# # plot b_c


# plot b_c
plt.figure()
for i in range(len(lw_list)):
    label = 'lw = '+str(lw_list[i])
    plt.errorbar(lc_list, y=b_c_data_avg[i,:], yerr = b_c_data_err[i,:], label=label)
plt.legend(fontsize=10)
plt.xlabel("l_c", fontsize=15)
plt.ylabel("b_c", fontsize=15)
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.tight_layout()
plt.savefig("b_c.PNG", dpi=300)
# plot b_c


# plot b_c (avg)
plt.figure()
max_vals = np.max(b_c_data_avg + b_c_data_err, axis=0)
min_vals = np.min(b_c_data_avg - b_c_data_err, axis=0)
err_vals = (max_vals - min_vals) / 2
# err_vals = np.mean(b_c_difference_err, axis=0)
plt.errorbar(lc_list, y=np.mean(b_c_data_avg, axis=0), yerr = err_vals, label='average', color='k',markersize=5, fmt='o')
log_x = np.log10(lc_list)
log_y = np.log10(np.mean(b_c_data_avg, axis=0))
log_sig = err_vals / (np.log(10) * np.mean(b_c_data_avg, axis=0))
w = 1 / log_sig**2
max_index = 3
try:
    coef = np.polyfit(log_x[0:max_index], log_y[0:max_index], deg=1, w=w[0:max_index])   # highest power first
    lc_fit_plot = np.linspace(lc_list[0], lc_list[max_index], 1000)
    bc_fit_plot = (lc_fit_plot**coef[0])*(10**coef[1])
    plt.plot(lc_fit_plot, bc_fit_plot, color='r', linestyle='--')
except:
    pass
plt.legend(fontsize=10)
plt.xlabel("l_c", fontsize=15)
plt.ylabel("b_c (avg)", fontsize=15)
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title("slope: "+str(round(coef[0],2)))
plt.tight_layout()
plt.savefig("b_c_avg.PNG", dpi=300)
# plot b_c (avg)



# plot b_c_difference
plt.figure()
for i in range(len(lc_list)):
    label = 'lw = '+str(lc_list[i])
    plt.errorbar(lc_list, y=b_c_difference_avg[i,:], yerr = b_c_difference_err[i,:], label=label)
plt.legend(fontsize=10)
plt.xlabel("l_c", fontsize=15)
plt.ylabel("b_c_difference", fontsize=15)
plt.yscale("log")
plt.xscale("log")
# plt.xlim((10, 120))
plt.grid()
plt.tight_layout()
plt.savefig("b_c_difference.PNG", dpi=300)
# plot b_c_difference


# plot b_c-b_c_terminal (avg)
plt.figure()
max_vals = np.max(b_c_difference_avg + b_c_difference_err, axis=0)
min_vals = np.min(b_c_difference_avg - b_c_difference_err, axis=0)
err_vals = (max_vals - min_vals) / 2
# err_vals = np.mean(b_c_difference_err, axis=0)
plt.errorbar(lc_list, y=np.mean(b_c_difference_avg, axis=0), yerr = err_vals, label='average', color='k',markersize=5, fmt='o')
log_x = np.log10(lc_list)
log_y = np.log10(np.mean(b_c_difference_avg, axis=0))
log_sig = err_vals / (np.log(10) * np.mean(b_c_difference_avg, axis=0))
w = 1 / log_sig**2
min_index = 2
max_index = 4
try:
    coef = np.polyfit(log_x[min_index:max_index], log_y[min_index:max_index], deg=1, w=w[min_index:max_index])   # highest power first
    lc_fit_plot = np.linspace(lc_list[min_index], lc_list[max_index], 1000)
    bc_diff_fit_plot = (lc_fit_plot**coef[0])*(10**coef[1])
    plt.plot(lc_fit_plot, bc_diff_fit_plot, color='r', linestyle='--')
except:
    pass
plt.legend(fontsize=10)
plt.xlabel("l_c", fontsize=15)
plt.ylabel("b_c_difference (avg)", fontsize=15)
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.title("slope: "+str(round(coef[0],2)))
plt.tight_layout()
plt.savefig("b_c_difference_avg.PNG", dpi=300)
# plot b_c-b_c_terminal (avg)




# plot b_c_ratio
plt.figure()
for i in range(len(lw_list)):
    label = 'lw = '+str(lw_list[i])
    plt.errorbar(lc_list, y=b_c_ratio_avg[i,:], yerr = b_c_ratio_err[i,:], label=label)
plt.legend(fontsize=10)
plt.xlabel("l_c", fontsize=15)
plt.ylabel("b_c_ratio", fontsize=15)
plt.yscale("log")
plt.xscale("log")
plt.grid()
plt.tight_layout()
plt.savefig("b_c_ratio.PNG", dpi=300)
# plot b_c_ratio

# # plot tau_w
# plt.figure()
# for i in range(len(lc_list)):
#     label = 'lc = '+str(lc_list[i])
#     plt.errorbar(lw_list, y=tau_w_data_avg[:,i], yerr = tau_w_data_err[:,i], label=label)
# plt.legend(fontsize=15)
# plt.xlabel("l_w", fontsize=15)
# plt.ylabel("tau_w", fontsize=15)
# plt.grid()
# plt.tight_layout()
# plt.savefig("tau_w.PNG", dpi=300)
# # plot tau_w

# # plot tau_c
# plt.figure()
# for i in range(len(lc_list)):
#     label = 'lc = '+str(lc_list[i])
#     plt.errorbar(lw_list, y=tau_c_data_avg[:,i], yerr = tau_c_data_err[:,i], label=label)
# plt.legend(fontsize=15)
# plt.xlabel("l_w", fontsize=15)
# plt.ylabel("tau_c", fontsize=15)
# plt.grid()
# plt.tight_layout()
# plt.savefig("tau_c.PNG", dpi=300)
# # plot tau_w

beta_aff_norm_plotter_both()


time_package.sleep(5)
script = "folder_zipper.sh"
title = "folder_zipper"
# Run script in a new terminal window
subprocess.run([
"gnome-terminal",
"--title=" + title,
"--", "bash", "-c",
f"./{script}; exec bash"
])
