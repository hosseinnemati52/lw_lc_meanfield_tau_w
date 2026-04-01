# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 11:55:30 2025

@author: nemat002
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
import os
import shutil
import subprocess

class paramsClass:
  # l_w_0 = 0 # how far the effect ON wt cells can reach
  # l_c_0 = 0 # how far the effect ON ca cells can reach
  test = 0

def reset_and_fill_folder(target_dir, sources):
    """
    target_dir: str
        Path to the folder you want to clear and repopulate.
    sources: list of paths
        Each path can be a file OR a folder to be copied into target_dir.
    """

    # 1. Remove the folder if it exists
    if os.path.exists(target_dir):
        # shutil.rmtree(target_dir)
        a =1

    # 2. Recreate it empty
    os.makedirs(target_dir, exist_ok=True)

    # 3. Copy each source (file or folder) into the target_dir
    for src in sources:
        if os.path.isdir(src):
            # copy folder → target_dir/src_name
            dst = os.path.join(target_dir, os.path.basename(src))
            shutil.copytree(src, dst)
        else:
            # copy file → target_dir/
            shutil.copy2(src, target_dir)

def save_object(obj, filename):
    """
    Save an object's attributes to a text file in JSON format.
    Only JSON-compatible attributes (numbers, strings, lists, dicts) will be saved.
    """
    with open(filename, "w") as f:
        json.dump(obj.__dict__, f, indent=4)

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

def params_updater(filepath, list_of_keys, list_of_vals):
    
    with open(filepath, "r") as f:
        data = json.load(f)

    # Modify values
    # data["l_w_0"] = lw*np.sqrt(data["a_w"])   # <-- your new value
    # data["l_c_0"] = lc*np.sqrt(data["a_c"])   # <-- your new value
    
    for i in range(len(list_of_keys)):
        key = list_of_keys[i]
        val = list_of_vals[i]
        data[key] = val   # <-- your new value
    
    # Save changes
    with open(filepath, "w") as f:
        json.dump(data, f, indent=4)
        
    return

def GD_init_vals_w():
    GD_init_intervals_w = np.loadtxt("GD_init_intervals_w.txt", delimiter='\t')
    
    if l_w >= GD_init_intervals_w[0,-1]:
        b_w_lower = GD_init_intervals_w[1,-1]
        b_w_upper = GD_init_intervals_w[2,-1]
        tau_w_lower = GD_init_intervals_w[3,-1]
        tau_w_upper = GD_init_intervals_w[4,-1]
    elif l_w <= GD_init_intervals_w[0,0]:
        b_w_lower = GD_init_intervals_w[1,0]
        b_w_upper = GD_init_intervals_w[2,0]
        tau_w_lower = GD_init_intervals_w[3,0]
        tau_w_upper = GD_init_intervals_w[4,0]
    else:
        b_w_lower = np.interp(l_w, GD_init_intervals_w[0,:], GD_init_intervals_w[1,:])
        b_w_upper = np.interp(l_w, GD_init_intervals_w[0,:], GD_init_intervals_w[2,:])
        tau_w_lower = np.interp(l_w, GD_init_intervals_w[0,:], GD_init_intervals_w[3,:])
        tau_w_upper = np.interp(l_w, GD_init_intervals_w[0,:], GD_init_intervals_w[4,:])
    
    b_w   = np.random.uniform(b_w_lower, b_w_upper)
    tau_w =  np.random.uniform(tau_w_lower, tau_w_upper)
    
    return b_w, tau_w

params = load_object(paramsClass, "params.txt")
l_w = params.l_w_0

# np.savetxt('GD_vals.csv', X=np.array([[l_w, b_w, b_c]]), delimiter=' , ', fmt='%.8e')
# initial random values

# initial random values
# b_w = np.random.uniform(0.005, 0.05)
# tau_w = np.random.uniform(1, 10)
b_w, tau_w = GD_init_vals_w()
b_c = np.random.uniform(1, 5)
#l_w = np.random.uniform(10, 100)
# tau_c = int(np.random.uniform(2000,5000))
# np.savetxt('GD_vals.csv', X=[b_w, b_c, tau_w, tau_c], delimiter=' , ', fmt='%.5e')

params_updater("params.txt", ["l_w_0", "b_w", "b_c", "tau_w"], [l_w, b_w, b_c, tau_w])

# initial run
subprocess.run(["python", "mean_field_lw_lc_ratio_multi.py"])  # run code here
# initial run

# reading the history
GD_LOG = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)
if GD_LOG.ndim == 1:
    # Convert to a 1-row matrix
    GD_LOG = GD_LOG.reshape(1, -1)
# b_w; b_c; tau_w; tau_c;
GD_vals = GD_LOG[-1,3:]
l_w = GD_vals[0]
b_w = GD_vals[1]
b_c = GD_vals[2]
tau_w = GD_vals[3]
# tau_c = GD_vals[3]
# reading the history


# make folder GD_temp
target = "GD_temp"
folders_to_copy = ["exp_data"]
files_to_copy = [
    "frame_switch.txt",
    "GD_log.csv",
    "mean_field_lw_lc_ratio_multi.py",
    "pp_plotter_multi.py", 
    "mixed_sample_bank.csv",
    "n_init_samples.csv",
    "params.txt",
    "PbyP_switch.txt", 
    "sample_indices_mix.csv",
    "y_d1_w_semi_exp.csv",
    "y_d2_w_semi_exp.csv",
    "y_d1_c_semi_exp.csv",
    "y_d2_c_semi_exp.csv",
    "lambdas.txt",
    "t0_k_coefs.csv"
]
sources = folders_to_copy + files_to_copy
reset_and_fill_folder(target, sources)
# make folder GD_temp


converge_cond = 0
n_iter_check = 20
conv_thresh = 0.01
time_thresh = 1 # 1 hour
sequence_check = dict()
sequence_check['cost'] = 1.0*np.ones(n_iter_check)
sequence_check['lw'] = 1.0*np.ones(n_iter_check)
sequence_check['bw'] = 1.0*np.ones(n_iter_check)
sequence_check['bc'] = 1.0*np.ones(n_iter_check)
sequence_check['tau_w'] = 1.0*np.ones(n_iter_check)
# sequence_check['tau_c'] = 1.0*np.ones(n_iter_check)

# GD variables
log_bw = np.log(b_w)
log_bc = np.log(b_c)
log_tau_w = np.log(tau_w)
# log_tau_c = np.log(tau_c)
# GD variables

# GD params
grad = dict()
grad['lw'] = 0.0
grad['log_bw'] = 0.0
grad['log_bc'] = 0.0
grad['log_tau_w'] = 0.0
# grad['log_tau_c'] = 0.0

lw_factor = 10.0
log_bw_factor = 1.0
log_bc_factor = 1.0
log_tau_w_factor = 3.0
# log_tau_c_factor = 1.0

w = 0.5  # update weight

pos_neg_perc = 0.005
delta_tau = 0.1 # 0.1 hour
late_learning_rate = 3e-3
init_learning_rate = 1.5e-3
# GD params

lw_list = [l_w]
bw_list = [b_w]
bc_list = [b_c]
tau_w_list = [tau_w]
# tau_c_list = [tau_c]

mother_dir = os.getcwd()

counter = 0
switch_bc_converge = 0 
# cost_estimate = 10**3
# cost_thresh_rate_change = 30*(l_w > 15.0) + 80*(l_w < 15.0)
while (not converge_cond):
    
    # learning_rate = late_learning_rate - (late_learning_rate-init_learning_rate)/(1+0.5*counter)
    
    learning_rate = init_learning_rate + (late_learning_rate-init_learning_rate)* ( 1.4*counter/(5 + 0.1 * counter**2) )
    
    #if counter<10:
    #    learning_rate = init_learning_rate
    #elif cost_estimate > cost_thresh_rate_change:
    #    learning_rate = 2 * late_learning_rate
    #else:
    #    learning_rate = 1 * late_learning_rate
    # learning_rate = 2e-3    
    ## grad evaluation
    
    # to l_w
    # #pos
    # l_w_pos = l_w + pos_neg_perc*abs(l_w)
    # params_updater("GD_temp/params.txt", ["l_w_0"], [l_w_pos])
    # os.chdir("GD_temp")   # go into daughter folder
    # subprocess.run(["python", "mean_field_lw_lc_ratio_multi.py"])  # run code here
    # try:
    #     cost_l_w_pos = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)[-1,0]
    # except IndexError:
    #     cost_l_w_pos = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)[0]
    # os.chdir(mother_dir)
    # #pos
    # #neg
    # l_w_neg = l_w - pos_neg_perc*abs(l_w)
    # params_updater("GD_temp/params.txt", ["l_w_0"], [l_w_neg])
    # os.chdir("GD_temp")   # go into daughter folder
    # subprocess.run(["python", "mean_field_lw_lc_ratio_multi.py"])  # run code here
    # try:
    #     cost_l_w_neg = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)[-1,0]
    # except IndexError:
    #     cost_l_w_neg = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)[0]
    # os.chdir(mother_dir)
    # #neg
    # grad['lw'] = (cost_l_w_pos-cost_l_w_neg)/(l_w_pos-l_w_neg)
    grad['lw'] = 0.0
    # to l_w
    
    # to log_bw
    #pos
    shutil.copy2("params.txt", "GD_temp/params.txt")
    log_bw_pos = log_bw + pos_neg_perc*abs(log_bw)
    b_w_pos = np.exp(log_bw_pos)
    # GD_vals_temp = np.array([[b_w_pos, b_c, tau_w, tau_c]])
    # GD_vals_temp = np.array([[b_w_pos, b_c]])
    # np.savetxt('GD_temp/'+'GD_vals.csv', X=GD_vals_temp, delimiter=' , ', fmt='%.8e')
    params_updater("GD_temp/params.txt", ["b_w"], [b_w_pos])
    os.chdir("GD_temp")   # go into daughter folder
    subprocess.run(["python", "mean_field_lw_lc_ratio_multi.py"])  # run code here
    try:
        cost_b_w_pos = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)[-1,0]
    except IndexError:
        cost_b_w_pos = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)[0]
    os.chdir(mother_dir)
    #pos
    #neg
    shutil.copy2("params.txt", "GD_temp/params.txt")
    log_bw_neg = log_bw - pos_neg_perc*abs(log_bw)
    b_w_neg = np.exp(log_bw_neg)
    # GD_vals_temp = np.array([[b_w_neg, b_c, tau_w, tau_c]])
    # GD_vals_temp = np.array([[b_w_neg, b_c]])
    # np.savetxt('GD_temp/'+'GD_vals.csv', X=GD_vals_temp, delimiter=' , ', fmt='%.8e')
    params_updater("GD_temp/params.txt", ["b_w"], [b_w_neg])
    os.chdir("GD_temp")   # go into daughter folder
    subprocess.run(["python", "mean_field_lw_lc_ratio_multi.py"])  # run code here
    try:
        cost_b_w_neg = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)[-1,0]
    except IndexError:
        cost_b_w_neg = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)[0]
    os.chdir(mother_dir)
    #neg
    grad['log_bw'] = (cost_b_w_pos-cost_b_w_neg)/(log_bw_pos-log_bw_neg)
    # to log_bw
    
    cost_estimate = (cost_b_w_pos+cost_b_w_neg)/2
    
    if switch_bc_converge:
        grad['log_bc'] = 0.0
    else:
        # to log_bc
        #pos
        shutil.copy2("params.txt", "GD_temp/params.txt")
        log_bc_pos = log_bc + pos_neg_perc*abs(log_bc)
        b_c_pos = np.exp(log_bc_pos)
        # GD_vals_temp = np.array([[b_w, b_c_pos, tau_w, tau_c]])
        # GD_vals_temp = np.array([[b_w, b_c_pos]])
        # np.savetxt('GD_temp/'+'GD_vals.csv', X=GD_vals_temp, delimiter=' , ', fmt='%.8e')
        params_updater("GD_temp/params.txt", ["b_c"], [b_c_pos])
        os.chdir("GD_temp")   # go into daughter folder
        subprocess.run(["python", "mean_field_lw_lc_ratio_multi.py"])  # run code here
        try:
            cost_b_c_pos = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)[-1,0]
        except IndexError:
            cost_b_c_pos = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)[0]
        os.chdir(mother_dir)
        #pos
        #neg
        shutil.copy2("params.txt", "GD_temp/params.txt")
        log_bc_neg = log_bc - pos_neg_perc*abs(log_bc)
        b_c_neg = np.exp(log_bc_neg)
        # GD_vals_temp = np.array([[b_w, b_c_neg, tau_w, tau_c]])
        # GD_vals_temp = np.array([[b_w, b_c_neg]])
        # np.savetxt('GD_temp/'+'GD_vals.csv', X=GD_vals_temp, delimiter=' , ', fmt='%.8e')
        params_updater("GD_temp/params.txt", ["b_c"], [b_c_neg])
        os.chdir("GD_temp")   # go into daughter folder
        subprocess.run(["python", "mean_field_lw_lc_ratio_multi.py"])  # run code here
        try:
            cost_b_c_neg = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)[-1,0]
        except IndexError:
            cost_b_c_neg = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)[0]
        os.chdir(mother_dir)
        #neg
        grad['log_bc'] = (cost_b_c_pos-cost_b_c_neg)/(log_bc_pos-log_bc_neg)
        # to log_bc
    
    # to log_tau_w
    #pos
    shutil.copy2("params.txt", "GD_temp/params.txt")
    # log_tau_w_pos = log_tau_w + pos_neg_perc*abs(log_tau_w)
    log_tau_w_pos = np.log(tau_w + delta_tau)
    tau_w_pos = np.exp(log_tau_w_pos)
    params_updater("GD_temp/params.txt", ["tau_w"], [tau_w_pos])
    os.chdir("GD_temp")   # go into daughter folder
    subprocess.run(["python", "mean_field_lw_lc_ratio_multi.py"])  # run code here
    try:
        cost_tau_w_pos = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)[-1,0]
    except IndexError:
        cost_tau_w_pos = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)[0]
    os.chdir(mother_dir)
    #pos
    #neg
    shutil.copy2("params.txt", "GD_temp/params.txt")
    # log_tau_w_neg = log_tau_w - pos_neg_perc*abs(log_tau_w)
    log_tau_w_neg = np.log(tau_w - delta_tau)
    tau_w_neg = np.exp(log_tau_w_neg)
    params_updater("GD_temp/params.txt", ["tau_w"], [tau_w_neg])
    os.chdir("GD_temp")   # go into daughter folder
    subprocess.run(["python", "mean_field_lw_lc_ratio_multi.py"])  # run code here
    try:
        cost_tau_w_neg = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)[-1,0]
    except IndexError:
        cost_tau_w_neg = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)[0]
    os.chdir(mother_dir)
    #neg
    grad['log_tau_w'] = (cost_tau_w_pos-cost_tau_w_neg)/(log_tau_w_pos-log_tau_w_neg)
    # to log_tau_w
    
    # to log_tau_c
    #pos
    # log_tau_c_pos = log_tau_c + pos_neg_perc*abs(log_tau_c)
    # tau_c_pos = np.exp(log_tau_c_pos)
    # GD_vals_temp = np.array([[b_w, b_c, tau_w, tau_c_pos]])
    # np.savetxt('GD_temp/'+'GD_vals.csv', X=GD_vals_temp, delimiter=' , ', fmt='%.5e')
    # os.chdir("GD_temp")   # go into daughter folder
    # subprocess.run(["python", "mean_field_lw_lc_ratio_multi.py"])  # run code here
    # try:
    #     cost_tau_c_pos = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)[-1,0]
    # except IndexError:
    #     cost_tau_c_pos = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)[0]
    # os.chdir(mother_dir)
    # #pos
    # #neg
    # log_tau_c_neg = log_tau_c - pos_neg_perc*abs(log_tau_c)
    # tau_c_neg = np.exp(log_tau_c_neg)
    # GD_vals_temp = np.array([[b_w, b_c, tau_w, tau_c_neg]])
    # np.savetxt('GD_temp/'+'GD_vals.csv', X=GD_vals_temp, delimiter=' , ', fmt='%.5e')
    # os.chdir("GD_temp")   # go into daughter folder
    # subprocess.run(["python", "mean_field_lw_lc_ratio_multi.py"])  # run code here
    # try:
    #     cost_tau_c_neg = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)[-1,0]
    # except IndexError:
    #     cost_tau_c_neg = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)[0]
    # os.chdir(mother_dir)
    # #neg
    # grad['log_tau_c'] = (cost_tau_c_pos-cost_tau_c_neg)/(log_tau_c_pos-log_tau_c_neg)
    # grad['log_tau_c'] = 0
    # to log_tau_c
    ## grad evaluation
        
    
    ## update variables
    delta_l_w =    - w * learning_rate * lw_factor * grad['lw']
    delta_log_bw = - w * learning_rate * log_bw_factor * grad['log_bw']
    delta_log_bc = - w * learning_rate * log_bc_factor * grad['log_bc']
    delta_log_tau_w = - w * learning_rate * log_tau_w_factor * grad['log_tau_w']
    # delta_log_tau_c = - w * learning_rate * log_tau_c_factor * grad['log_tau_c']
    
    l_w += delta_l_w
    log_bw += delta_log_bw
    log_bc += delta_log_bc
    log_tau_w += delta_log_tau_w
    # log_tau_c += delta_log_tau_c
    
    b_w  = np.exp(log_bw)
    b_c  = np.exp(log_bc)
    tau_w  = np.exp(log_tau_w)
    # tau_c  = np.exp(log_tau_c)
    
    
    # GD_vals = np.array([[b_w, b_c, tau_w, tau_c]])
    # GD_vals = np.array([[b_w, b_c]])
    # np.savetxt('GD_vals.csv', X=GD_vals, delimiter=' , ', fmt='%.8e')
    
    # params_updater("params.txt", ["l_w_0", "b_w", "b_c"], [l_w, b_w, b_c])
    params_updater("params.txt", ["l_w_0", "b_w", "b_c", "tau_w"], [l_w, b_w, b_c, tau_w])
    # new run
    subprocess.run(["python", "mean_field_lw_lc_ratio_multi.py"])  # run code here
    # new run
    new_cost = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)[-1,0]
    
    lw_list.append(l_w)
    bw_list.append(b_w)
    bc_list.append(b_c)
    tau_w_list.append(tau_w)
    # tau_c_list.append(tau_c)
    ## update variables

    
    ## updating sequences
    sequence_check['cost'][:-1] = sequence_check['cost'][1:]
    sequence_check['cost'][-1]  = new_cost # new
    sequence_check['lw'][:-1] = sequence_check['lw'][1:]
    sequence_check['lw'][-1]  = l_w # new
    sequence_check['bw'][:-1] = sequence_check['bw'][1:]
    sequence_check['bw'][-1]  = b_w # new
    sequence_check['bc'][:-1] = sequence_check['bc'][1:]
    sequence_check['bc'][-1]  = b_c # new
    sequence_check['tau_w'][:-1] = sequence_check['tau_w'][1:]
    sequence_check['tau_w'][-1]  = tau_w # new
    # sequence_check['tau_c'][:-1] = sequence_check['tau_c'][1:]
    # sequence_check['tau_c'][-1]  = tau_c # new
    ## updating sequences
    
    ## converge condition
    cost_mat_check = np.abs(sequence_check['cost'][:,None]-sequence_check['cost'][None,:]) \
        /np.abs(sequence_check['cost'][:,None])
    lw_mat_check = np.abs(sequence_check['lw'][:,None]-sequence_check['lw'][None,:]) \
        /np.abs(sequence_check['lw'][:,None])
    bw_mat_check = np.abs(sequence_check['bw'][:,None]-sequence_check['bw'][None,:]) \
        /np.abs(sequence_check['bw'][:,None])
    bc_mat_check = np.abs(sequence_check['bc'][:,None]-sequence_check['bc'][None,:]) \
        /np.abs(sequence_check['bc'][:,None])
    # tau_w_mat_check = np.abs(sequence_check['tau_w'][:,None]-sequence_check['tau_w'][None,:]) \
    #     /np.abs(sequence_check['tau_w'][:,None])
    tau_w_mat_check = np.abs(sequence_check['tau_w'][:,None]-sequence_check['tau_w'][None,:])
    
    # converge_cond = \
    #     np.max(lw_mat_check) < conv_thresh and \
    #     np.max(bw_mat_check) < conv_thresh and \
    #     np.max(bc_mat_check) < conv_thresh and \
    #     np.max(tau_w_mat_check) < time_thresh
    #     # np.max(tau_w_mat_check) < conv_thresh
    converge_cond = (np.max(cost_mat_check) < conv_thresh)
        
        
        
        
    if np.max(bc_mat_check) < conv_thresh:
        switch_bc_converge = 1
    # converge_cond = \
    # max(abs(np.diff(sequence_check['bw'])/sequence_check['bw'][:-1])) < conv_thresh and \
    # max(abs(np.diff(sequence_check['bc'])/sequence_check['bc'][:-1])) < conv_thresh #and \
    # max(abs(np.diff(sequence_check['tau_w'])/sequence_check['tau_w'][:-1])) < conv_thresh and \
    # max(abs(np.diff(sequence_check['tau_c'])/sequence_check['tau_c'][:-1])) < conv_thresh
    ## converge condition
    
    # # new run
    # subprocess.run(["python", "mean_field_lw_lc_ratio_multi.py"])  # run code here
    # # new run
    
    counter += 1
    
    print("********************************")
    print("iteration: "+str(counter))
    print("********************************")
    
print("##################")
print("CONVERGENCE: pp running!")
subprocess.run(["python", "cost_mat_calc.py"])
subprocess.run(["python", "pp_plotter_multi.py"])
print("##################")
    





