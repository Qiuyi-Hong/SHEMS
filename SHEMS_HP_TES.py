"""
This file contains the code for modelling SHEMS integrated with HP and TES. 
"""

import pyomo.environ as pyo 

model = pyo.AbstractModel(name="SHEMS_HP_TES")

############################################################################
######################## Setting parameters: ###############################
############################################################################
model.T = pyo.RangeSet(1, 24, 1)

# Electricity demand of the house in kW
model.d_ele = pyo.Param(model.T, mutable=True)

# Heat demand of the house in kW
model.d_heat = pyo.Param(model.T, mutable=True)

# Electricity import price in p/kWh
model.pi_import = pyo.Param(model.T, mutable=True)

# Electricity export price in p/kWh
model.pi_export = pyo.Param(model.T, mutable=True)

# Battery energy storage params:
# model.c_BESS = pyo.Param()
model.eta_c_BESS = pyo.Param()
model.eta_d_BESS = pyo.Param()
model.epsilon_BESS = pyo.Param() # Self-discharge rate
model.E_min_BESS = pyo.Param()
model.E_max_BESS = pyo.Param()
model.E_init = pyo.Param()
model.p_c_min_BESS = pyo.Param() 
model.p_c_max_BESS = pyo.Param()
model.p_d_min_BESS = pyo.Param()
model.p_d_max_BESS = pyo.Param()
model.delta_t = pyo.Param()


# Air-to-water heat pump params:
model.q_HP_max = pyo.Param()
model.rho_HP = pyo.Param()
# model.p_HP = pyo.Param(model.T, mutable=True)
model.COP = pyo.Param(model.T, mutable=True)
# model.a = pyo.Param()
# model.b = pyo.Param()
# model.T_water = pyo.Param(model.T, mutable=True)
# model.T_air = pyo.Param(model.T, mutable=True)

# Thermal energy storage params:
model.eta_c_TES = pyo.Param()
model.eta_d_TES = pyo.Param()
model.epsilon_TES = pyo.Param()
model.Q_min_TES = pyo.Param()
model.Q_max_TES = pyo.Param()
model.Q_init_TES = pyo.Param()
model.q_c_min_TES = pyo.Param()
model.q_c_max_TES = pyo.Param()
model.q_d_min_TES = pyo.Param()
model.q_d_max_TES = pyo.Param()


# PV generation in kW
model.p_pv = pyo.Param(model.T, mutable=True)

############################################################################
######################## Setting decision variables: ####################### 
############################################################################

# Battery energy storage vars:
model.p_c_BESS = pyo.Var(model.T)
model.p_d_BESS = pyo.Var(model.T)
model.E_BESS = pyo.Var(model.T)
model.gamma_c = pyo.Var(model.T, within=pyo.Boolean)
model.gamma_d = pyo.Var(model.T, within=pyo.Boolean)

# Air-to-water heat pump vars:
model.q_HP = pyo.Var(model.T)
model.p_HP = pyo.Var(model.T)
model.sigma = pyo.Var(model.T, within=pyo.Boolean)

# Thermal energy storage vars:
model.q_c_TES = pyo.Var(model.T)
model.q_d_TES = pyo.Var(model.T)
model.Q_TES = pyo.Var(model.T)
model.theta_c = pyo.Var(model.T, within=pyo.Boolean)
model.theta_d = pyo.Var(model.T, within=pyo.Boolean)

# Electricity volume imported from the grid in kWh
model.p_import = pyo.Var(model.T, within=pyo.NonNegativeReals)

# Electricity volume exported to the grid in kWh
model.p_export = pyo.Var(model.T, within=pyo.NonNegativeReals)

############################################################################
######################## Setting constraints: ############################## 
############################################################################

# Battery energy storage constraints:
def energyConstr(model, t):
    if t == 1:
        return model.E_BESS[t] == model.E_init + model.eta_c_BESS * model.p_c_BESS[t] * model.delta_t - (1/model.eta_d_BESS) * model.p_d_BESS[t] * model.delta_t - model.epsilon_BESS * model.delta_t
    else:
        return model.E_BESS[t] == model.E_BESS[t-1] + model.eta_c_BESS * model.p_c_BESS[t] * model.delta_t - (1/model.eta_d_BESS) * model.p_d_BESS[t] * model.delta_t - model.epsilon_BESS * model.delta_t

model.energyConstr = pyo.Constraint(model.T, rule=energyConstr)

def energyMax(model, t):
    return model.E_BESS[t] <= model.E_max_BESS 

model.energyMax = pyo.Constraint(model.T, rule=energyMax)

def energyMin(model, t):
    return model.E_BESS[t] >= model.E_min_BESS 

model.energyMin = pyo.Constraint(model.T, rule=energyMin)

def energyEquivalence(model):
    return model.E_init == model.E_BESS[24]

model.energyEquivalence = pyo.Constraint(rule=energyEquivalence)

def powerChargeMax(model, t):
    return model.p_c_BESS[t] <= model.gamma_c[t] * model.p_c_max_BESS 

model.powerChargeMax = pyo.Constraint(model.T, rule=powerChargeMax)

def powerChargeMin(model, t):
    return model.p_c_BESS[t] >= model.gamma_c[t] * model.p_c_min_BESS 

model.powerChargeMin = pyo.Constraint(model.T, rule=powerChargeMin)

def powerDischargeMax(model, t):
    return model.p_d_BESS[t] <= model.gamma_d[t] * model.p_d_max_BESS 

model.powerDischargeMax = pyo.Constraint(model.T, rule=powerDischargeMax)

def powerDischargeMin(model, t):
    return model.p_d_BESS[t] >= model.gamma_d[t] * model.p_d_min_BESS 

model.powerDischargeMin = pyo.Constraint(model.T, rule=powerDischargeMin)

def gammaConstr(model, t):
    return model.gamma_c[t] + model.gamma_d[t] <= 1 

model.gammaConstr = pyo.Constraint(model.T, rule=gammaConstr)

# Air-to-water heat pump constraints:
# def COPVal(model, t):
#     return model.COP[t] == model.a * (model.T_water[t] - model.T_air[t]) + model.b

# model.COPVal = pyo.Constraint(model.T, rule=COPVal)

def COPDefinition(model, t):
    return model.p_HP[t] * model.COP[t] == model.q_HP[t]

model.COPDefinition = pyo.Constraint(model.T, rule=COPDefinition)

def powerHPMin(model, t):
    return model.rho_HP * model.q_HP_max * model.sigma[t] <= model.q_HP[t]

model.powerHPMin = pyo.Constraint(model.T, rule=powerHPMin)

def powerHPMax(model, t):
    return model.q_HP[t] <= model.q_HP_max * model.sigma[t]

model.powerHPMax = pyo.Constraint(model.T, rule=powerHPMax)

# Thermal energy storage constraints:
def energyConstrTES(model, t):
    if t == 1:
        return model.Q_TES[t] == model.Q_init_TES + model.eta_c_TES * model.q_c_TES[t] * model.delta_t - (1/model.eta_d_TES) * model.q_d_TES[t] * model.delta_t - model.epsilon_TES * model.delta_t
    else:
        return model.Q_TES[t] == model.Q_TES[t-1] + model.eta_c_TES * model.q_c_TES[t] * model.delta_t - (1/model.eta_d_TES) * model.q_d_TES[t] * model.delta_t - model.epsilon_TES * model.delta_t
    
model.energyConstrTES = pyo.Constraint(model.T, rule=energyConstrTES)

def chargingeqv(model, t):
    return model.q_c_TES[t] <= model.q_HP[t]

model.chargingeqv = pyo.Constraint(model.T, rule=chargingeqv)

def energyTESBalance(model):
    return model.Q_init_TES == model.Q_TES[24]

model.energyTESBalance = pyo.Constraint(rule=energyTESBalance)

def energyMinTES(model, t):
    return model.Q_min_TES <= model.Q_TES[t]

model.energyMinTES = pyo.Constraint(model.T, rule=energyMinTES)

def energyMaxTES(model, t):
    return model.Q_TES[t] <= model.Q_max_TES

model.energyMaxTES = pyo.Constraint(model.T, rule=energyMaxTES)

def powerChargeMinTES(model, t):
    return model.theta_c[t] * model.q_c_min_TES <= model.q_c_TES[t]

model.powerChargeMinTES = pyo.Constraint(model.T, rule=powerChargeMinTES)

def powerChargeMaxTES(model, t):
    return model.q_c_TES[t] <= model.theta_c[t] * model.q_c_max_TES

model.powerChargeMaxTES = pyo.Constraint(model.T, rule=powerChargeMaxTES)

def powerDischargeMinTES(model, t):
    return model.theta_d[t] * model.q_d_min_TES <= model.q_d_TES[t]

model.powerDischargeMinTES = pyo.Constraint(model.T, rule=powerDischargeMinTES)

def powerDischargeMaxTES(model, t):
    return model.q_d_TES[t] <= model.theta_d[t] * model.q_d_max_TES

model.powerDischargeMaxTES = pyo.Constraint(model.T, rule=powerDischargeMaxTES)

def thetaConstr(model, t):
    return model.theta_c[t] + model.theta_d[t] <= 1

model.thetaConstr = pyo.Constraint(model.T, rule=thetaConstr)

# Power balance constraints:
def powerBalance(model, t):
    return model.p_HP[t] * model.delta_t + model.d_ele[t] + model.p_c_BESS[t] * model.delta_t + model.p_export[t] == model.p_d_BESS[t] * model.delta_t + model.p_pv[t] * model.delta_t + model.p_import[t]

model.powerBalance = pyo.Constraint(model.T, rule=powerBalance)

# Heat balance constraints:
def heatBalance(model, t):
    return model.d_heat[t] == model.q_d_TES[t] * model.delta_t

model.heatBalance = pyo.Constraint(model.T, rule=heatBalance)

############################################################################
######################## Setting objective function: ####################### 
############################################################################

def ObjectiveFuction(model):
    total = 0 
    for t in model.T:
        total += model.pi_import[t] * model.p_import[t] - model.pi_export[t] * model.p_export[t]
    return total

model.obj = pyo.Objective(rule=ObjectiveFuction, sense=pyo.minimize)