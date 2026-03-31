#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 17 23:57:34 2025

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


# exp data
overal_WT_mix = np.loadtxt("exp_data/"+"overal_WT_mix.csv", delimiter=',')
overal_C_mix = np.loadtxt("exp_data/"+"overal_C_mix.csv", delimiter=',')
exp_time   = overal_WT_mix[0,:]
WT_mix_norm_avg = overal_WT_mix[1,:]
WT_mix_norm_err = overal_WT_mix[2,:]
C_mix_norm_avg  = overal_C_mix[1,:]
C_mix_norm_err  = overal_C_mix[2,:]
# exp data

W_mix_fit_smples_coefs = np.loadtxt("W_mix_fit_smples_coefs.csv", delimiter=',')
C_mix_fit_smples_coefs = np.loadtxt("C_mix_fit_smples_coefs.csv", delimiter=',')

# WT
a_list = W_mix_fit_smples_coefs[0,:]
b_list = W_mix_fit_smples_coefs[1,:]
c_list = W_mix_fit_smples_coefs[2,:]

random_sample_ind_w = np.random.randint(np.shape(W_mix_fit_smples_coefs)[1])
target_a = a_list[random_sample_ind_w]
target_b = b_list[random_sample_ind_w]
target_c = c_list[random_sample_ind_w]

d1_mat = np.zeros((len(a_list), len(exp_time)))
d2_mat = np.zeros((len(a_list), len(exp_time)))
for i in range(len(a_list)):
    d1_mat[i,:] = 2*a_list[i]*exp_time+b_list[i]
    d2_mat[i,:] = 2*a_list[i]
    
# discrete_y_d1_avg_exp = np.mean(d1_mat, axis=0)
discrete_y_d1_avg_exp = 2*target_a*exp_time+target_b
# discrete_y_d1_err_exp = np.std(d1_mat, axis=0)/np.sqrt(len(a_list))
discrete_y_d1_err_exp = np.std(d1_mat, axis=0)
# discrete_y_d2_avg_exp = np.mean(d2_mat, axis=0)
discrete_y_d2_avg_exp = 2 * target_a
# discrete_y_d2_err_exp = np.std(d2_mat, axis=0)/np.sqrt(len(a_list))
discrete_y_d2_err_exp = np.std(d2_mat, axis=0)
target_abc_fit = np.array([random_sample_ind_w, target_a, target_b, target_c])
np.savetxt("target_abc_fit.csv", target_abc_fit, fmt='%.6f', delimiter=',')



y_d1_w_semi_exp = np.zeros((3,len(exp_time)))
y_d2_w_semi_exp = np.zeros((3,len(exp_time)))

y_d1_w_semi_exp[0,:] = exp_time
y_d1_w_semi_exp[1,:] = discrete_y_d1_avg_exp
y_d1_w_semi_exp[2,:] = discrete_y_d1_err_exp

y_d2_w_semi_exp[0,:] = exp_time
y_d2_w_semi_exp[1,:] = discrete_y_d2_avg_exp
y_d2_w_semi_exp[2,:] = discrete_y_d2_err_exp

np.savetxt("y_d1_w_semi_exp.csv", y_d1_w_semi_exp, fmt='%.5e', delimiter=',')
np.savetxt("y_d2_w_semi_exp.csv", y_d2_w_semi_exp, fmt='%.5e', delimiter=',')
# WT


pure_beta_mat = np.loadtxt("pure_organoid_exponents.csv", delimiter=',')
beta_w = pure_beta_mat[0,0]
beta_w_err = pure_beta_mat[1,0]
beta_c = pure_beta_mat[0,1]
beta_c_err = pure_beta_mat[1,1]

# C
L_list  = C_mix_fit_smples_coefs[0,:]
t0_list = C_mix_fit_smples_coefs[1,:]
k_list  = C_mix_fit_smples_coefs[2,:]
y0_list = C_mix_fit_smples_coefs[3,:]

random_sample_ind_c = np.random.randint(np.shape(C_mix_fit_smples_coefs)[1])
target_L  =  L_list[random_sample_ind_c]
target_t0 = t0_list[random_sample_ind_c]
target_k  =  k_list[random_sample_ind_c]
target_y0 = y0_list[random_sample_ind_c]
target_coefs_c_fit = np.array([random_sample_ind_c, target_L, target_t0, target_k, target_y0])
np.savetxt("target_coefs_c_fit.csv", target_coefs_c_fit, fmt='%.6f', delimiter=',')


d1_mat = np.zeros((len(L_list), len(exp_time)))
d2_mat = np.zeros((len(L_list), len(exp_time)))
for i in range(len(L_list)):
    expo = np.exp(-k_list[i]*(exp_time-t0_list[i]))
    d1_mat[i,:] = k_list[i]*L_list[i]*expo / (1+expo)**2
    d2_mat[i,:] = k_list[i]*k_list[i]*L_list[i]*expo * (-1+expo) / (1+expo)**3
d1_mat = d1_mat + beta_c

# discrete_y_d1_avg_exp = np.mean(d1_mat, axis=0)
expo = np.exp(-target_k*(exp_time-target_t0))
discrete_y_d1_avg_exp = beta_c + target_k*target_L*expo / (1+expo)**2
# discrete_y_d1_err_exp = np.std(d1_mat, axis=0)/np.sqrt(len(L_list))
discrete_y_d1_err_exp = np.std(d1_mat, axis=0)
discrete_y_d1_err_exp = (discrete_y_d1_err_exp**2 + beta_c_err**2)**0.5

# discrete_y_d2_avg_exp = np.mean(d2_mat, axis=0)
expo = np.exp(-target_k*(exp_time-target_t0))
discrete_y_d2_avg_exp = target_k*target_k*target_L*expo * (-1+expo) / (1+expo)**3
# discrete_y_d2_err_exp = np.std(d2_mat, axis=0)/np.sqrt(len(L_list))
discrete_y_d2_err_exp = np.std(d2_mat, axis=0)

y_d1_c_semi_exp = np.zeros((3,len(exp_time)))
y_d2_c_semi_exp = np.zeros((3,len(exp_time)))

y_d1_c_semi_exp[0,:] = exp_time
y_d1_c_semi_exp[1,:] = discrete_y_d1_avg_exp
y_d1_c_semi_exp[2,:] = discrete_y_d1_err_exp

y_d2_c_semi_exp[0,:] = exp_time
y_d2_c_semi_exp[1,:] = discrete_y_d2_avg_exp
y_d2_c_semi_exp[2,:] = discrete_y_d2_err_exp

np.savetxt("y_d1_c_semi_exp.csv", y_d1_c_semi_exp, fmt='%.5e', delimiter=',')
np.savetxt("y_d2_c_semi_exp.csv", y_d2_c_semi_exp, fmt='%.5e', delimiter=',')
# C