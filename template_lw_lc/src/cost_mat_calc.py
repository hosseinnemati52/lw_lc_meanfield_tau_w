# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 16:38:47 2025

@author: Nemat002
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

class paramsClass:
  # l_w_0 = 0 # how far the effect ON wt cells can reach
  # l_c_0 = 0 # how far the effect ON ca cells can reach
  test = 0
  
class initClass:
  n_w_init = 0
  n_c_init = 0

def params_func():
    params = paramsClass()
    params.char_cell_size = 10 #micron 
    params.l_w_0 = 1.0 * params.char_cell_size # how far the effect ON wt cells can reach
    params.l_c_0 = 1.0 * params.char_cell_size # how far the effect ON ca cells can reach
    
    params.beta_w_unaff = 0.0284
    params.beta_c_unaff = 0.0398
    
    params.beta_w_aff = 0.01
    params.beta_c_aff = 0.0398
    
    return params

def init_func(params, init_n):
    
    init = initClass()
    
    init.n_w_init = init_n.n_w_init
    init.n_c_init = init_n.n_c_init
    
    init.A_w_init = init.n_w_init * params.a_w
    init.A_c_init = init.n_c_init * params.a_c
    
    total_A = init.A_w_init + init.A_c_init
    
    init.r = np.sqrt(total_A/(4*np.pi))
    
    # A = 2*\pi*(r^2)*(1 - cos \theta)
    init.theta_2 = np.arccos( 1 - init.A_c_init / (2 * np.pi * init.r**2) ) # border of wt and c
    
    # split of c
    init.theta_1 = max(0, init.theta_2 - params.l_c_0/init.r)
    init.A_c_unaff = 2 * np.pi * (init.r**2) * (1-np.cos(init.theta_1))
    init.A_c_aff   = init.A_c_init - init.A_c_unaff
    # if init.theta_1>0:
    #     init.A_c_unaff = 2 * np.pi * (init.r**2) * (1-np.cos(init.theta_1))
    #     init.A_c_aff   = init.A_c_init - init.A_c_unaff
    # elif init.theta_1 <= 0:
    #     init.theta_1 = 0
    #     init.A_c_unaff = 0
    #     init.A_c_aff   = init.A_c_init
    # split of c

    # split of wt    
    init.theta_3 = min(np.pi, init.theta_2 + params.l_w_0/init.r)
    init.A_w_unaff = 2 * np.pi * (init.r**2) * (1-np.cos(np.pi - init.theta_3))
    init.A_w_aff   = init.A_w_init - init.A_w_unaff
    # if init.theta_3 < np.pi:
    #     init.A_w_unaff = 2 * np.pi * (init.r**2) * (1-np.cos(np.pi - init.theta_3))
    #     init.A_w_aff   = init.A_w_init - init.A_w_unaff
    # elif init.theta_3 >= np.pi:
    #     init.theta_3 = np.pi
    #     init.A_w_unaff = 0
    #     init.A_w_aff   = init.A_w_init
    # split of wt
    
    return init


def init_numbers_maker():
    
    try:
        n_init_samples = np.loadtxt('n_init_samples.csv', dtype=int, delimiter=',')
        n_org = np.shape(n_init_samples)[0]
    except:
        n_org = 200
        mixed_sample_bank = np.loadtxt("mixed_sample_bank.csv", delimiter=',', dtype=int)
        size = np.shape(mixed_sample_bank)[0]
        sample_indices_mix = np.random.randint(0,size,n_org)
        np.savetxt('sample_indices_mix.csv', X=sample_indices_mix, delimiter=',', fmt='%d')
        
        n_init_samples = np.zeros((n_org,2), dtype=int)
        for org_c in range(n_org):
            n_init_samples[org_c,0] = mixed_sample_bank[sample_indices_mix[org_c],0]
            n_init_samples[org_c,1] = mixed_sample_bank[sample_indices_mix[org_c],1]
        
        np.savetxt('n_init_samples.csv', X=n_init_samples, fmt='%d', delimiter=',')
    
    return n_init_samples

def cost_calc(cost_key):
    
    if cost_key == "w":
        factor_w = 1
        factor_c = 0
    elif cost_key == "c":
        factor_w = 0
        factor_c = 1
    elif cost_key == "both":
        factor_w = 1
        factor_c = 1
        
        
    cost = 0.0
    
    log_w_bar_mod = np.log(np.mean(A_w_mat / A_w_mat[:, [0]], axis=0))
    log_w_bar_mod_err =  (np.std(A_w_mat / A_w_mat[:, [0]], axis=0)/np.sqrt(n_org)) / (np.mean(A_w_mat / A_w_mat[:, [0]], axis=0))
    
    log_w_bar_exp     = np.log(WT_mix_norm_avg)
    log_w_bar_exp_err = WT_mix_norm_err / WT_mix_norm_avg
    
    log_c_bar_mod = np.log(np.mean(A_c_mat / A_c_mat[:, [0]], axis=0))
    log_c_bar_mod_err =  (np.std(A_c_mat / A_c_mat[:, [0]], axis=0)/np.sqrt(n_org)) / (np.mean(A_c_mat / A_c_mat[:, [0]], axis=0))
    
    log_c_bar_exp = np.log(C_mix_norm_avg)
    log_c_bar_exp_err = C_mix_norm_err / C_mix_norm_avg
    
    # plt.figure()
    # plt.errorbar(exp_time, y=log_w_bar_exp, yerr = log_w_bar_exp_err, color='m')
    # plt.errorbar(exp_time, y=log_c_bar_exp, yerr = log_c_bar_exp_err, color='g')
    # plt.errorbar(time, y=log_w_bar_mod, yerr = log_w_bar_mod_err, color='r')
    # plt.errorbar(time, y=log_c_bar_mod, yerr = log_c_bar_mod_err, color='b')
    
    
    log_w_bar_mod_interpolate     = np.interp(exp_time, time, log_w_bar_mod)
    log_w_bar_mod_interpolate_err = np.interp(exp_time, time, log_w_bar_mod_err)
    
    log_c_bar_mod_interpolate     = np.interp(exp_time, time, log_c_bar_mod)
    log_c_bar_mod_interpolate_err = np.interp(exp_time, time, log_c_bar_mod_err)
    
    cost_weights_w = 1/log_w_bar_exp_err**2
    cost_weights_c = 1/log_c_bar_exp_err**2
    
    cost_w = float(np.sum(cost_weights_w[1:] * (log_w_bar_mod_interpolate[1:]-log_w_bar_exp[1:])**2))
    cost_c = float(np.sum(cost_weights_c[1:] * (log_c_bar_mod_interpolate[1:]-log_c_bar_exp[1:])**2))
    
    cost = dict()
    
    cost['w'] = factor_w * cost_w
    cost['c'] = factor_c * cost_c
    cost['tot'] = cost_w + cost_c
    
    return cost

def cost_calc_derivs(cost_key):
    
    if cost_key == "w":
        factor_w = 1
        factor_c = 0
    elif cost_key == "c":
        factor_w = 0
        factor_c = 1
    elif cost_key == "both":
        factor_w = 1
        factor_c = 1
        
    lambdas = np.loadtxt("lambdas.txt", delimiter=',')
    cost_mat = np.zeros((2,3))
    # each row for a type (WT or C)
    # each col for y or y' or y''
    # cost = lambda_0 \sum cost_y + lambda_1 \sum cost_y' + lambda_2 \sum cost_y''
    
    n_points_exp = len(WT_mix_norm_avg)
    ############## WT ###################
    # W_mix_fit_smples_coefs = np.loadtxt("W_mix_fit_smples_coefs.csv", delimiter=',')
    # C_mix_fit_smples_coefs = np.loadtxt("C_mix_fit_smples_coefs.csv", delimiter=',')
    
    # y
    log_w_bar_mod = np.log(np.mean(A_w_mat / A_w_mat[:, [0]], axis=0))
    log_w_bar_mod_err =  (np.std(A_w_mat / A_w_mat[:, [0]], axis=0)/np.sqrt(n_org)) / (np.mean(A_w_mat / A_w_mat[:, [0]], axis=0))
    log_w_bar_mod_interpolate     = np.interp(exp_time, time, log_w_bar_mod)
    log_w_bar_mod_interpolate_err = np.interp(exp_time, time, log_w_bar_mod_err)
    log_w_bar_exp     = np.log(WT_mix_norm_avg)
    log_w_bar_exp_err = WT_mix_norm_err / WT_mix_norm_avg
    # SEM_y = (  log_w_bar_mod_interpolate_err**2 + log_w_bar_exp_err**2  )**0.5
    SEM_y = log_w_bar_exp_err
    r_y = abs(log_w_bar_mod_interpolate-log_w_bar_exp)
    # cost_mat[0,0] = lambdas[0] * np.sum( (r_y[1:]/SEM_y[1:])**2 )
    cost_mat[0,0] = np.sum( (r_y[1:]/SEM_y[1:])**2 )
    # y
    
    # y' and y''
    spline_w = UnivariateSpline(time, log_w_bar_mod, s=0)
    deriv_1_w = spline_w.derivative(1)(time)
    deriv_2_w = spline_w.derivative(2)(time)
    discrete_y_d1_avg = np.interp(exp_time, time, deriv_1_w) # avg for y' at the descrete points: t=0, 5, 10, ..., 70
    discrete_y_d2_avg = np.interp(exp_time, time, deriv_2_w) # avg for y'' at the descrete points: t=0, 5, 10, ..., 70
    n_blocks = 10
    n_org_per_block = int(n_org/n_blocks)
    discrete_y_d1_mat = np.zeros((n_blocks, len(exp_time)))
    discrete_y_d2_mat = np.zeros((n_blocks, len(exp_time)))
    for block_c in range(n_blocks):
        beg_ind = block_c * n_org_per_block
        end_ind = beg_ind + n_org_per_block
        A_w_mat_block = A_w_mat[beg_ind:end_ind,:]
        log_w_bar_mod_block = np.log(np.mean(A_w_mat_block / A_w_mat_block[:, [0]], axis=0))
        spline_w_block = UnivariateSpline(time, log_w_bar_mod_block, s=0)
        deriv_1_w_block = spline_w_block.derivative(1)(time)
        deriv_2_w_block = spline_w_block.derivative(2)(time)
        discrete_y_d1_mat[block_c,:] = np.interp(exp_time, time, deriv_1_w_block)
        discrete_y_d2_mat[block_c,:] = np.interp(exp_time, time, deriv_2_w_block)
    discrete_y_d1_err = np.std(discrete_y_d1_mat, axis=0)/np.sqrt(n_blocks-1) # SEM for y' at the descrete points: t=0, 5, 10, ..., 70
    discrete_y_d2_err = np.std(discrete_y_d2_mat, axis=0)/np.sqrt(n_blocks-1) # SEM for y'' at the descrete points: t=0, 5, 10, ..., 70
    
    # a_list = W_mix_fit_smples_coefs[0,:]
    # b_list = W_mix_fit_smples_coefs[1,:]
    # c_list = W_mix_fit_smples_coefs[2,:]
    # d1_mat = np.zeros((len(a_list), len(exp_time)))
    # d2_mat = np.zeros((len(a_list), len(exp_time)))
    # for i in range(len(a_list)):
    #     d1_mat[i,:] = 2*a_list[i]*exp_time+b_list[i]
    #     d2_mat[i,:] = 2*a_list[i]
    # discrete_y_d1_avg_exp = np.mean(d1_mat, axis=0)
    # discrete_y_d1_err_exp = np.std(d1_mat, axis=0)/np.sqrt(len(a_list))
    # discrete_y_d2_avg_exp = np.mean(d2_mat, axis=0)
    # discrete_y_d2_err_exp = np.std(d2_mat, axis=0)/np.sqrt(len(a_list))
    try:
        y_d1_w_semi_exp = np.loadtxt("y_d1_w_semi_exp.csv", delimiter=',')
        y_d2_w_semi_exp = np.loadtxt("y_d2_w_semi_exp.csv", delimiter=',')
    except FileNotFoundError:
        subprocess.run(["python", "semi_exp_deriv.py"])
        y_d1_w_semi_exp = np.loadtxt("y_d1_w_semi_exp.csv", delimiter=',')
        y_d2_w_semi_exp = np.loadtxt("y_d2_w_semi_exp.csv", delimiter=',')
    discrete_y_d1_avg_exp = y_d1_w_semi_exp[1,:]
    discrete_y_d1_err_exp = y_d1_w_semi_exp[2,:]
    discrete_y_d2_avg_exp = y_d2_w_semi_exp[1,:]
    discrete_y_d2_err_exp = y_d2_w_semi_exp[2,:]
    
    SEM_y_d1 = (  discrete_y_d1_err**2 + discrete_y_d1_err_exp**2  )**0.5
    r_y_d1 = abs(discrete_y_d1_avg-discrete_y_d1_avg_exp)
    # cost_mat[0,1] = lambdas[1] * np.sum( (r_y_d1[1:]/SEM_y_d1[1:])**2 )
    cost_mat[0,1] = np.sum( (r_y_d1[1:]/SEM_y_d1[1:])**2 )
    
    SEM_y_d2 = (  discrete_y_d2_err**2 + discrete_y_d2_err_exp**2  )**0.5
    r_y_d2 = abs(discrete_y_d2_avg-discrete_y_d2_avg_exp)
    # cost_mat[0,2] = lambdas[2] * np.sum( (r_y_d2[1:]/SEM_y_d2[1:])**2 )
    cost_mat[0,2] = np.sum( (r_y_d2[1:]/SEM_y_d2[1:])**2 )
    # y' and y''
    ############## WT ###################
    
    ############## C ###################
    # y
    log_c_bar_mod = np.log(np.mean(A_c_mat / A_c_mat[:, [0]], axis=0))
    log_c_bar_mod_err =  (np.std(A_c_mat / A_c_mat[:, [0]], axis=0)/np.sqrt(n_org)) / (np.mean(A_c_mat / A_c_mat[:, [0]], axis=0))
    log_c_bar_mod_interpolate     = np.interp(exp_time, time, log_c_bar_mod)
    log_c_bar_mod_interpolate_err = np.interp(exp_time, time, log_c_bar_mod_err)
    log_c_bar_exp = np.log(C_mix_norm_avg)
    log_c_bar_exp_err = C_mix_norm_err / C_mix_norm_avg
    # SEM_y = (  log_c_bar_mod_interpolate_err**2 + log_c_bar_exp_err**2  )**0.5
    SEM_y = log_c_bar_exp_err
    r_y = abs(log_c_bar_mod_interpolate-log_c_bar_exp)
    # cost_mat[1,0] = lambdas[0] * np.sum( (r_y[1:]/SEM_y[1:])**2 )
    cost_mat[1,0] = np.sum( (r_y[1:]/SEM_y[1:])**2 )
    # y
    
    # y' and y''
    spline_c = UnivariateSpline(time, log_c_bar_mod, s=0)
    deriv_1_c = spline_c.derivative(1)(time)
    deriv_2_c = spline_c.derivative(2)(time)
    discrete_y_d1_avg = np.interp(exp_time, time, deriv_1_c) # avg for y' at the descrete points: t=0, 5, 10, ..., 70
    discrete_y_d2_avg = np.interp(exp_time, time, deriv_2_c) # avg for y'' at the descrete points: t=0, 5, 10, ..., 70
    n_blocks = 10
    n_org_per_block = int(n_org/n_blocks)
    discrete_y_d1_mat = np.zeros((n_blocks, len(exp_time)))
    discrete_y_d2_mat = np.zeros((n_blocks, len(exp_time)))
    for block_c in range(n_blocks):
        beg_ind = block_c * n_org_per_block
        end_ind = beg_ind + n_org_per_block
        A_c_mat_block = A_c_mat[beg_ind:end_ind,:]
        log_c_bar_mod_block = np.log(np.mean(A_c_mat_block / A_c_mat_block[:, [0]], axis=0))
        spline_c_block = UnivariateSpline(time, log_c_bar_mod_block, s=0)
        deriv_1_c_block = spline_c_block.derivative(1)(time)
        deriv_2_c_block = spline_c_block.derivative(2)(time)
        discrete_y_d1_mat[block_c,:] = np.interp(exp_time, time, deriv_1_c_block)
        discrete_y_d2_mat[block_c,:] = np.interp(exp_time, time, deriv_2_c_block)
    discrete_y_d1_err = np.std(discrete_y_d1_mat, axis=0)/np.sqrt(n_blocks-1) # SEM for y' at the descrete points: t=0, 5, 10, ..., 70
    discrete_y_d2_err = np.std(discrete_y_d2_mat, axis=0)/np.sqrt(n_blocks-1) # SEM for y'' at the descrete points: t=0, 5, 10, ..., 70
    
    y_d1_c_semi_exp = np.loadtxt("y_d1_c_semi_exp.csv", delimiter=',')
    y_d2_c_semi_exp = np.loadtxt("y_d2_c_semi_exp.csv", delimiter=',')
    discrete_y_d1_avg_exp = y_d1_c_semi_exp[1,:]
    discrete_y_d1_err_exp = y_d1_c_semi_exp[2,:]
    discrete_y_d2_avg_exp = y_d2_c_semi_exp[1,:]
    discrete_y_d2_err_exp = y_d2_c_semi_exp[2,:]
    
    SEM_y_d1 = (  discrete_y_d1_err**2 + discrete_y_d1_err_exp**2  )**0.5
    r_y_d1 = abs(discrete_y_d1_avg-discrete_y_d1_avg_exp)
    # cost_mat[1,1] = lambdas[1] * np.sum( (r_y_d1[1:]/SEM_y_d1[1:])**2 )
    cost_mat[1,1] = np.sum( (r_y_d1[1:]/SEM_y_d1[1:])**2 )
    
    SEM_y_d2 = (  discrete_y_d2_err**2 + discrete_y_d2_err_exp**2  )**0.5
    r_y_d2 = abs(discrete_y_d2_avg-discrete_y_d2_avg_exp)
    # cost_mat[1,2] = lambdas[2] * np.sum( (r_y_d2[1:]/SEM_y_d2[1:])**2 )
    cost_mat[1,2] = np.sum( (r_y_d2[1:]/SEM_y_d2[1:])**2 )
    # y' and y''
    ############## C ###################
    
    # cost_w = np.sum(cost_mat[0,:])
    # cost_c = np.sum(cost_mat[1,:])
    
    cost_w = np.sum(cost_mat[0,:] * lambdas)
    cost_c = np.sum(cost_mat[1,:] * lambdas)
    
    cost = dict()
    
    cost['w'] = factor_w * cost_w
    cost['c'] = factor_c * cost_c
    cost['tot'] = cost_w + cost_c
    
    return cost, cost_mat

def GD_logger():
    
    try:
        # cost_tot; cost_w; cost_c; b_w; b_c; tau_w; tau_c;
        GD_log = np.loadtxt("GD_log.csv", delimiter=',', dtype=float)
        if GD_log.ndim == 1:
            # Convert to a 1-row matrix
            GD_log = GD_log.reshape(1, -1)
            
        GD_log_new = np.vstack([GD_log, np.zeros((1, GD_log.shape[1]))])
        # GD_log_new[-1,:] = np.array([cost_dict['tot'], cost_dict['w'], cost_dict['c'], b_w, b_c, tau_w, tau_c])
        # GD_log_new[-1,:] = np.array([cost_dict['tot'], cost_dict['w'], cost_dict['c'], l_w, b_w, b_c])
        GD_log_new[-1,:] = np.array([cost_dict['tot'], cost_dict['w'], cost_dict['c'], l_w, b_w, b_c, tau_w])
        np.savetxt('GD_log.csv', X=GD_log_new, delimiter=' , ', fmt='%.6e')
    except:
        GD_log_new = np.zeros((1,7))
        # GD_log_new = np.zeros((1,6))
        # GD_log_new[-1,:] = np.array([cost_dict['tot'], cost_dict['w'], cost_dict['c'], b_w, b_c, tau_w, tau_c])
        # GD_log_new[-1,:] = np.array([cost_dict['tot'], cost_dict['w'], cost_dict['c'], l_w, b_w, b_c])
        GD_log_new[-1,:] = np.array([cost_dict['tot'], cost_dict['w'], cost_dict['c'], l_w, b_w, b_c, tau_w])
        np.savetxt('GD_log.csv', X=GD_log_new, delimiter=' , ', fmt='%.6e')
    
    return 0

def cost_logger(cost_matrix):
    
    with open("cost_log.txt", "a") as f:
        for row in cost_matrix:
            f.write(" ".join(map(str, row)) + "\n")
        f.write("--------\n")
    
    return 0
    
# exp data
overal_WT_mix = np.loadtxt("exp_data/"+"overal_WT_mix.csv", delimiter=',')
overal_C_mix = np.loadtxt("exp_data/"+"overal_C_mix.csv", delimiter=',')
exp_time   = overal_WT_mix[0,:]
WT_mix_norm_avg = overal_WT_mix[1,:]
WT_mix_norm_err = overal_WT_mix[2,:]
C_mix_norm_avg  = overal_C_mix[1,:]
C_mix_norm_err  = overal_C_mix[2,:]
# exp data

# sim params
time = np.loadtxt("data/"+'time.txt',  delimiter=',')
# sim params


# reading params
# params = params_func()
# save_object(params, "params.txt")
params = load_object(paramsClass, "params.txt")
# reading params

# temporal delta functions
# b_w = 0.013; tau_w = 20;
# # b_c = 0.1; tau_c = 80;
# # b_w = 0.3; tau_w = 20;
# b_c = 2.0; tau_c = 70;

# GD_vals = np.loadtxt("GD_vals.csv", delimiter=',')
# b_w; b_c; tau_w; tau_c;
# b_w = GD_vals[0]
# b_c = GD_vals[1]
b_w = params.b_w
b_c = params.b_c
beta_w_term = params.beta_w_term_ratio * params.beta_w_unaff
tau_w = params.tau_w
# tau_w = GD_vals[2]
# tau_c = GD_vals[3]

# delta_w = b_w * (1-np.exp(-time/tau_w))
# delta_c = b_c * (1-np.exp(-time/tau_c))

# delta_w = b_w * time
# delta_c = b_c * time
t0_k_coefs = np.loadtxt("t0_k_coefs.csv", delimiter = ',')
t0     = t0_k_coefs[1]
k_coef = t0_k_coefs[2]
# temporal delta functions

r_mat = np.loadtxt("data/"+'r_mat.txt', delimiter=',')
# np.savetxt("data/"+'th_1_mat.txt', th_1_mat, fmt='%.5e', delimiter=',')
# np.savetxt("data/"+'th_2_mat.txt', th_2_mat, fmt='%.5e', delimiter=',')
# np.savetxt("data/"+'th_3_mat.txt', th_3_mat, fmt='%.5e', delimiter=',')
# np.savetxt("data/"+'th_v_w_mat.txt', th_v_w_mat, fmt='%.5e', delimiter=',')
# np.savetxt("data/"+'th_v_c_mat.txt', th_v_c_mat, fmt='%.5e', delimiter=',')
A_w_mat = np.loadtxt("data/"+'A_w_mat.txt', delimiter=',')
A_c_mat = np.loadtxt("data/"+'A_c_mat.txt', delimiter=',')
n_org = np.shape(A_w_mat)[0]

# np.savetxt("data/"+'beta_w_aff_mat.txt', beta_w_aff_mat, fmt='%.5e', delimiter=',')
# np.savetxt("data/"+'beta_c_aff_mat.txt', beta_c_aff_mat, fmt='%.5e', delimiter=',')

# np.savetxt("data/"+'A_w_aff_mat.txt',   A_w_aff_mat,   fmt='%.5e', delimiter=',')
# np.savetxt("data/"+'A_w_unaff_mat.txt', A_w_unaff_mat, fmt='%.5e', delimiter=',')
# np.savetxt("data/"+'A_c_aff_mat.txt',   A_c_aff_mat,   fmt='%.5e', delimiter=',')
# np.savetxt("data/"+'A_c_unaff_mat.txt', A_c_unaff_mat, fmt='%.5e', delimiter=',')
# np.savetxt("data/"+'A_c_v_mat.txt',   A_c_v_mat,   fmt='%.5e', delimiter=',')
# np.savetxt("data/"+'A_w_v_mat.txt',   A_w_v_mat, fmt='%.5e', delimiter=',')

# cost_dict = cost_calc("both")
cost_dict, cost_mat = cost_calc_derivs("both")

np.savetxt("data/"+'cost_mat_final.txt', cost_mat, fmt='%.5e', delimiter=',')
