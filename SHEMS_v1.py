"""
This file contains the code for modelling SHEMS. 
"""

import pyomo.environ as pyo 

model = pyo.AbstractModel(name="SHEMS")

############################################################################
######################## Setting parameters: ###############################
############################################################################
model.T = pyo.RangeSet(1, 48, 1)
model.TT = pyo.RangeSet(0, 48, 1)

# Non-schedulable electricity demand of the house in kW
model.d_non = pyo.Param(model.T, mutable=True)

# electricity import price in p/kWh
model.pi_import = pyo.Param(model.T, mutable=True)

# electricity export price in p/kWh
model.pi_export = pyo.Param(model.T, mutable=True)

# Energy storage params:
model.c_ESS = pyo.Param()
model.eta_c_ESS = pyo.Param()
model.eta_d_ESS = pyo.Param()
model.epsilon_ESS = pyo.Param() # Self-discharge rate
model.E_min_ESS = pyo.Param()
model.E_max_ESS = pyo.Param()
model.E_init = pyo.Param()
model.p_c_min_ESS = pyo.Param() 
model.p_c_max_ESS = pyo.Param()
model.p_d_min_ESS = pyo.Param()
model.p_d_max_ESS = pyo.Param()
model.delta_t = pyo.Param()

# Schedulable appliances -- uninterruptible:
model.num_appliances = pyo.Param()
model.A = pyo.RangeSet(1, model.num_appliances)
model.N_a = pyo.Param(model.A, within=pyo.NonNegativeReals)
model.t_start = pyo.Param(model.A, within=pyo.NonNegativeReals)
model.t_stop = pyo.Param(model.A, within=pyo.NonNegativeReals)
model.p_app = pyo.Param(model.A, within=pyo.NonNegativeReals)

# Self-Owned EV:
model.c_EV = pyo.Param()
model.eta_c_EV = pyo.Param()
model.eta_d_EV = pyo.Param()
model.epsilon_EV = pyo.Param() # Self-discharge rate
model.E_min_EV = pyo.Param()
model.E_max_EV = pyo.Param()
model.p_c_min_EV = pyo.Param() 
model.p_c_max_EV = pyo.Param()
model.p_d_min_EV = pyo.Param()
model.p_d_max_EV = pyo.Param()
model.E_arr_EV = pyo.Param()

model.t_arr_EV = pyo.Param()
model.t_dep_EV = pyo.Param()

# External EV:
model.num_external_EV = pyo.Param()
model.N = pyo.RangeSet(1, model.num_external_EV)
model.eta_c_EEV = pyo.Param(model.N, within=pyo.NonNegativeReals)
model.eta_d_EEV = pyo.Param(model.N, within=pyo.NonNegativeReals)
model.epsilon_EEV = pyo.Param(model.N, within=pyo.NonNegativeReals)
model.alpha_EEV = pyo.Param(model.N, within=pyo.NonNegativeReals)
model.E_max_EEV = pyo.Param(model.N, within=pyo.NonNegativeReals)
model.E_min_EEV = pyo.Param(model.N, within=pyo.NonNegativeReals)
model.E_arr_EEV = pyo.Param(model.N, within=pyo.NonNegativeReals)
model.p_c_min_EEV = pyo.Param(model.N, within=pyo.NonNegativeReals)
model.p_c_max_EEV = pyo.Param(model.N, within=pyo.NonNegativeReals)

model.t_arr_EEV = pyo.Param(model.N, within=pyo.NonNegativeReals)
model.t_dep_EEV = pyo.Param(model.N, within=pyo.NonNegativeReals)

# PV generation in kW
model.q_pv = pyo.Param(model.T, mutable=True)

############################################################################
######################## Setting decision variables: ####################### 
############################################################################

# Energy storage vars:
model.p_c_ESS = pyo.Var(model.T)
model.p_d_ESS = pyo.Var(model.T)
model.E_ESS = pyo.Var(model.T)
model.gamma_c = pyo.Var(model.T, within=pyo.Boolean)
model.gamma_d = pyo.Var(model.T, within=pyo.Boolean)

# Schedulable appliances vars:
model.x = pyo.Var(model.TT, model.A, within=pyo.Boolean)
model.phi = pyo.Var(model.TT, model.A, within=pyo.Boolean)

# Self-Owned EV vars:
model.p_c_EV = pyo.Var(model.T)
model.p_d_EV = pyo.Var(model.T)
model.E_EV = pyo.Var(model.T)
model.zeta_c = pyo.Var(model.T, within=pyo.Boolean)
model.zeta_d = pyo.Var(model.T, within=pyo.Boolean)

# Self-Owned EV vars:
model.p_c_EEV = pyo.Var(model.T, model.N)
model.p_d_EEV = pyo.Var(model.T, model.N)
model.E_EEV = pyo.Var(model.T, model.N)

# electricity volume imported from the grid in kWh
model.q_import = pyo.Var(model.T, within=pyo.NonNegativeReals)

# electricity volume exported to the grid in kWh
model.q_export = pyo.Var(model.T, within=pyo.NonNegativeReals)

############################################################################
######################## Setting constraints: ############################## 
############################################################################

# Energy storage constraints:
def energyConstr(model, t):
    if t == 1:
        return model.E_ESS[t] == model.E_init + model.eta_c_ESS * model.p_c_ESS[t] * model.delta_t - (1/model.eta_d_ESS) * model.p_d_ESS[t] * model.delta_t - model.epsilon_ESS * model.delta_t
    else:
        return model.E_ESS[t] == model.E_ESS[t-1] + model.eta_c_ESS * model.p_c_ESS[t] * model.delta_t - (1/model.eta_d_ESS) * model.p_d_ESS[t] * model.delta_t - model.epsilon_ESS * model.delta_t

model.energyConstr = pyo.Constraint(model.T, rule=energyConstr)

def energyMax(model, t):
    return model.E_ESS[t] <= model.E_max_ESS 

model.energyMax = pyo.Constraint(model.T, rule=energyMax)

def energyMin(model, t):
    return model.E_ESS[t] >= model.E_min_ESS 

model.energyMin = pyo.Constraint(model.T, rule=energyMin)

def energyEquivalence(model):
    return model.E_init == model.E_ESS[48]

model.energyEquivalence = pyo.Constraint(rule=energyEquivalence)

def powerChargeMax(model, t):
    return model.p_c_ESS[t] <= model.gamma_c[t] * model.p_c_max_ESS 

model.powerChargeMax = pyo.Constraint(model.T, rule=powerChargeMax)

def powerChargeMin(model, t):
    return model.p_c_ESS[t] >= model.gamma_c[t] * model.p_c_min_ESS 

model.powerChargeMin = pyo.Constraint(model.T, rule=powerChargeMin)

def powerDischargeMax(model, t):
    return model.p_d_ESS[t] <= model.gamma_d[t] * model.p_d_max_ESS 

model.powerDischargeMax = pyo.Constraint(model.T, rule=powerDischargeMax)

def powerDischargeMin(model, t):
    return model.p_d_ESS[t] >= model.gamma_d[t] * model.p_d_min_ESS 

model.powerDischargeMin = pyo.Constraint(model.T, rule=powerDischargeMin)

def gammaConstr(model, t):
    return model.gamma_c[t] + model.gamma_d[t] <= 1 

model.gammaConstr = pyo.Constraint(model.T, rule=gammaConstr)

# Scheduleable appliances constraints:
def numConstr(model, a):
    return sum(model.x[t, a] for t in range(model.t_start[a], model.t_stop[a])) == model.N_a[a]

model.numConstr = pyo.Constraint(model.A, rule=numConstr)

def schedulableConstr(model, a):
    return sum(model.x[t, a] for t in range(1, model.t_start[a])) + sum(model.x[t, a] for t in range(model.t_stop[a], 48+1)) == 0

model.scheduleConstr = pyo.Constraint(model.A, rule=schedulableConstr)

def uninteruptibleConstr1(model, t, a):
    return model.x[t, a] <= 1 - model.phi[t, a]

model.uniteruptibleConstr1 = pyo.Constraint(model.T, model.A, rule=uninteruptibleConstr1)

def uninteruptibleConstr2(model, t, a):
    return model.x[t - 1, a] - model.x[t, a] <= model.phi[t, a] 

model.uniteruptibleConstr2 = pyo.Constraint(model.T, model.A, rule=uninteruptibleConstr2)

def uninteruptibleConstr3(model, t, a):
    return model.phi[t - 1, a] <= model.phi[t, a]

model.uninteruptibleConstr3 = pyo.Constraint(model.T, model.A, rule=uninteruptibleConstr3)

# Self-Owned EV constraints:
def EVEnergyConstr(model, t):
    if t == model.t_arr_EV:
        return model.E_EV[t] == model.E_arr_EV + model.eta_c_EV * model.p_c_EV[t] * model.delta_t - (1/model.eta_d_EV) * model.p_d_EV[t] * model.delta_t - model.epsilon_EV * model.delta_t
    else:
        return model.E_EV[t] == model.E_EV[t-1] + model.eta_c_EV * model.p_c_EV[t] * model.delta_t - (1/model.eta_d_EV) * model.p_d_EV[t] * model.delta_t - model.epsilon_EV * model.delta_t
    
model.EVEnergyConstr = pyo.Constraint(pyo.RangeSet(model.t_arr_EV, model.t_dep_EV - 1), rule=EVEnergyConstr)

def EVDepConstr(model):
    return model.E_max_EV * 0.8 <= model.E_EV[model.t_dep_EV]

model.EVDepConstr = pyo.Constraint(rule=EVDepConstr)

def EVMaxEnergyConstr(model, t):
    return model.E_EV[t] <= model.E_max_EV

model.EVMaxEnergyConstr = pyo.Constraint(pyo.RangeSet(model.t_arr_EV, model.t_dep_EV), rule=EVMaxEnergyConstr)

def EVMinEnergyConstr(model, t):
    return model.E_EV[t] >= model.E_min_EV

model.EVMinEnergyConstr = pyo.Constraint(pyo.RangeSet(model.t_arr_EV, model.t_dep_EV), rule=EVMinEnergyConstr)

def EVChargeMax(model, t):
    return model.p_c_EV[t] <= model.zeta_c[t] * model.p_c_max_EV

model.EVChargeMax = pyo.Constraint(pyo.RangeSet(model.t_arr_EV, model.t_dep_EV - 1), rule=EVChargeMax)

def EVChargeMin(model, t):
    return model.p_c_EV[t] >= model.zeta_c[t] * model.p_c_min_EV

model.EVChargeMin = pyo.Constraint(pyo.RangeSet(model.t_arr_EV, model.t_dep_EV - 1), rule=EVChargeMin)

def EVDischargeMax(model, t):
    return model.p_d_EV[t] <= model.zeta_d[t] * model.p_d_max_EV

model.EVDischargeMax = pyo.Constraint(pyo.RangeSet(model.t_arr_EV, model.t_dep_EV - 1), rule=EVDischargeMax)

def EVDischargeMin(model, t):
    return model.p_d_EV[t] >= model.zeta_d[t] * model.p_d_min_EV

model.EVDischargeMin = pyo.Constraint(pyo.RangeSet(model.t_arr_EV, model.t_dep_EV - 1), rule=EVDischargeMin)

def EVChargeDischargeConstr(model):
    return sum(model.p_c_EV[t] + model.p_d_EV[t] for t in range(1, int(model.t_arr_EV))) + sum(model.p_c_EV[t] + model.p_d_EV[t] for t in range(int(model.t_dep_EV), 48+1)) == 0

model.EVChargeDischargeConstr = pyo.Constraint(rule=EVChargeDischargeConstr)

# External EV constraints:
t_arr_EEVEnergyConstr = 0
t_dep_EEVEnergyConstr = 0
def EEVEnergyConstr(model, t, n):
    global t_arr_EEVEnergyConstr, t_dep_EEVEnergyConstr
    t_arr_EEVEnergyConstr = model.t_arr_EEV[n]
    t_dep_EEVEnergyConstr = model.t_dep_EEV[n]
    if t == model.t_arr_EEV[n]:
        return model.E_EEV[t, n] == model.E_arr_EEV[n] + model.eta_c_EEV[n] * model.p_c_EEV[t, n] * model.delta_t - model.epsilon_EEV[n] * model.delta_t
    else:
        return model.E_EEV[t, n] == model.E_EEV[t-1, n] + model.eta_c_EEV[n] * model.p_c_EEV[t, n] * model.delta_t - model.epsilon_EEV[n] * model.delta_t
    
model.EEVEnergyConstr = pyo.Constraint(pyo.RangeSet(t_arr_EEVEnergyConstr, t_dep_EEVEnergyConstr - 1), model.N, rule=EEVEnergyConstr)

def EEVDepConstr(model, n):
    return model.E_max_EEV[n] * model.alpha[n] == model.E_EEV[model.t_dep_EEV, n]

model.EEVDepConstr = pyo.Constraint(model.N, rule=EEVDepConstr)

t_arr_EEVChargeMax = 0
t_dep_EEVChargeMax = 0
def EEVChargeMax(model, t, n):
    global t_arr_EEVChargeMax, t_dep_EEVChargeMax
    t_arr_EEVChargeMax = model.t_arr_EEV[n]
    t_dep_EEVChargeMax = model.t_dep_EEV[n]
    return model.p_c_EEV[t, n] <= model.p_c_max_EEV[n]

model.EEVChargeMax = pyo.Constraint(pyo.RangeSet(t_arr_EEVChargeMax, t_dep_EEVChargeMax - 1), model.N, rule=EEVChargeMax)

t_arr_EEVChargeMin = 0
t_dep_EEVChargeMin = 0
def EEVChargeMin(model, t, n):
    global t_arr_EEVChargeMin, t_dep_EEVChargeMin
    t_arr_EEVChargeMin = model.t_arr_EEV[n]
    t_dep_EEVChargeMin = model.t_dep_EEV[n]
    return model.p_c_EEV[t, n] >= model.p_c_min_EEV[n]

model.EEVChargeMin = pyo.Constraint(pyo.RangeSet(t_arr_EEVChargeMin, t_dep_EEVChargeMin - 1), model.N, rule=EEVChargeMin)

def EEVPriceConstr(model, t):
    return model.pi_EEV[t] == model.pi_export[t] * 1.2

model.EEVPriceConstr = pyo.Constraint(model.T, rule=EEVPriceConstr)

def EEVChargeDischargeConstr(model, n):
    return sum(model.p_c_EEV[t, n] for t in range(1, model.t_arr_EEV)) + sum(model.p_c_EEV[t, n] for t in range(model.t_dep_EEV, 48+1)) == 0

model.EEVChargeDischargeConstr = pyo.Constraint(model.N, rule=EEVChargeDischargeConstr)

# Power balance constraints:
def powerBalance(model, t):
    return sum(model.p_c_EEV[t,n] * model.delta_t for n in range(1, model.num_external_EV+1)) + model.p_c_EV[t] * model.delta_t + sum(model.p_app[a] * model.x[a, t] * model.delta_t for a in range(1, model.num_appliances+1)) + model.d_non[t] * model.delta_t + model.p_c_ESS[t] * model.delta_t + model.q_export[t] == model.p_d_ESS[t] * model.delta_t + model.q_pv[t] * model.delta_t + model.q_import[t] + model.p_d_EV[t] * model.delta_t


model.powerBalance = pyo.Constraint(model.T, rule=powerBalance)

############################################################################
######################## Setting objective function: ####################### 
############################################################################

def ObjectiveFuction(model):
    total = 0 
    for t in model.T:
        total += model.pi_import[t] * model.q_import[t] - model.pi_export[t] * model.q_export[t] + model.c_ESS * (model.p_c_ESS[t] + model.p_d_ESS[t]) * model.delta_t + model.c_EV * (model.p_c_EV[t] + model.p_d_EV[t]) * model.delta_t - sum(model.pi_EEV[t] * model.p_c_EEV[t, n] * model.delta_t for n in range(1, model.num_external_EV+1))
    return total

model.obj = pyo.Objective(rule=ObjectiveFuction, sense=pyo.minimize)