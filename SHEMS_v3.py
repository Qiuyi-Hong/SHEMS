"""
This file contains the code for modelling SHEMS v3. 
"""

import pyomo.environ as pyo 

model = pyo.AbstractModel(name="SHEMS_v3")

############################################################################
######################## Setting parameters: ###############################
############################################################################
model.t_end = pyo.Param(initialize=48)
model.T = pyo.RangeSet(1, model.t_end, 1)
model.delta_t = pyo.Param()

# Electricity demand of the house in kW
model.d_ele = pyo.Param(model.T, mutable=True)

# Heat demand of the house in kW
model.d_heat = pyo.Param(model.T, mutable=True)

# Electricity import price in p/kWh
model.pi_import = pyo.Param(model.T, mutable=True, within=pyo.Reals)

# Electricity export price in p/kWh
model.pi_export = pyo.Param(model.T, mutable=True, within=pyo.Reals)


# Air-to-water heat pump params:
model.T_out = pyo.Param(model.T, mutable=True)
model.COP = pyo.Param(model.T, mutable=True)
model.q_HP_max = pyo.Param()
model.q_HP_min = pyo.Param()

# Comfort params:
# For SH:
model.m_in = pyo.Param()
model.c_in = pyo.Param()
model.T_in_LB = pyo.Param()
model.T_in_UB = pyo.Param()
model.K_SH = pyo.Param()
model.T_in_init = pyo.Param()
# For DHW:
model.m_TES = pyo.Param()
model.c_TES = pyo.Param()
model.T_TES_LB = pyo.Param()
model.T_TES_UB = pyo.Param()
model.K_TES = pyo.Param()
model.T_TES_init = pyo.Param()

# Thermal energy storage params:
model.Q_TES_min = pyo.Param()
model.Q_TES_max = pyo.Param()
model.Q_TES_init = pyo.Param()
model.q_TES_d_min = pyo.Param()
model.q_TES_d_max = pyo.Param()
model.V_TES = pyo.Param()
model.rho_TES = pyo.Param()


# PV generation in kW
model.p_pv = pyo.Param(model.T, mutable=True)

############################################################################
######################## Setting decision variables: ####################### 
############################################################################

# Air-to-water heat pump vars:
model.T_TES = pyo.Var(model.T)
model.q_HP = pyo.Var(model.T)
model.p_HP = pyo.Var(model.T)
model.sigma_SH = pyo.Var(model.T, within=pyo.Boolean)
model.sigma_DHW = pyo.Var(model.T, within=pyo.Boolean)
model.q_SH = pyo.Var(model.T)
model.q_DHW = pyo.Var(model.T)

# Comfort vars:
# For SH:
model.T_in = pyo.Var(model.T)
model.epsilon_SH = pyo.Var(model.T)
# For DHW:
model.epsilon_TES = pyo.Var(model.T)

# Thermal energy storage vars:
model.Q_TES = pyo.Var(model.T)
model.q_TES_d = pyo.Var(model.T)
model.T_init = pyo.Var(model.T)


# Electricity volume imported from the grid in kWh
model.p_import = pyo.Var(model.T, within=pyo.NonNegativeReals)

# Electricity volume exported to the grid in kWh
model.p_export = pyo.Var(model.T, within=pyo.NonNegativeReals)

############################################################################
######################## Setting constraints: ############################## 
############################################################################

# Air-to-water heat pump constraints:
def COPDefinition(model, t):
    return model.p_HP[t] * model.COP[t] == model.q_HP[t]

model.COPDefinition = pyo.Constraint(model.T, rule=COPDefinition)

def powerSHMin(model, t):
    return model.sigma_SH[t] * model.q_HP_min <= model.q_SH[t]

model.powerSHMin = pyo.Constraint(model.T, rule=powerSHMin)

def powerSHMax(model, t):
    return model.q_SH[t] <= model.sigma_SH[t] * model.q_HP_max

model.powerSHMax = pyo.Constraint(model.T, rule=powerSHMax)

def powerDHWMin(model, t):
    return model.sigma_DHW[t] * model.q_HP_min <= model.q_DHW[t]

model.powerDHWMin = pyo.Constraint(model.T, rule=powerDHWMin)

def powerDHWMax(model, t):
    return model.q_DHW[t] <= model.sigma_DHW[t] * model.q_HP_max

model.powerDHWMax = pyo.Constraint(model.T, rule=powerDHWMax)

def sigmaConstr(model, t):
    return model.sigma_SH[t] + model.sigma_DHW[t] <= 1

model.sigmaConstr = pyo.Constraint(model.T, rule=sigmaConstr)

def SHDHWConstr(model, t):
    return model.q_SH[t] + model.q_DHW[t] == model.q_HP[t]

model.SHDHWConstr = pyo.Constraint(model.T, rule=SHDHWConstr)

# Comfort constraints:
# For SH:
def TempSH(model, t):
    if t == 1:
        return model.T_in[t] == model.T_in_init + ((model.q_SH[t] - model.epsilon_SH[t]) * model.delta_t) * 3.6e6 / (model.m_in * model.c_in)
    else:
        return model.T_in[t] == model.T_in[t-1] + ((model.q_SH[t] - model.epsilon_SH[t]) * model.delta_t) * 3.6e6 / (model.m_in * model.c_in)

model.TempSH = pyo.Constraint(model.T, rule=TempSH)

def TempconstrSHMin(model, t):
    return model.T_in_LB <= model.T_in[t]

model.TempconstrSHMin = pyo.Constraint(model.T, rule=TempconstrSHMin)

def TempconstrSHMax(model, t):
    return model.T_in[t] <= model.T_in_UB

model.TempconstrSHMax = pyo.Constraint(model.T, rule=TempconstrSHMax)

def lossSH(model, t):
    return model.epsilon_SH[t] == model.K_SH * (model.T_in[t] - model.T_out[t])

model.lossSH = pyo.Constraint(model.T, rule=lossSH)

def TempEqvSH(model):
    return model.T_in[model.t_end] == model.T_in_init

model.TempEqvSH = pyo.Constraint(rule=TempEqvSH)

# For DHW:
def TempDHW(model, t):
    if t == 1:
        return model.T_TES[t] == model.T_TES_init + ((model.q_DHW[t] - model.q_TES_d[t] - model.epsilon_TES[t]) * model.delta_t) * 3.6e6 / (model.m_TES * model.c_TES)
    else:
        return model.T_TES[t] == model.T_TES[t-1] + ((model.q_DHW[t] - model.q_TES_d[t] - model.epsilon_TES[t]) * model.delta_t) * 3.6e6 / (model.m_TES * model.c_TES)

model.TempDHW = pyo.Constraint(model.T, rule=TempDHW)

def TempconstrDHWMin(model, t):
    return model.T_TES_LB <= model.T_TES[t]

model.TempconstrDHWMin = pyo.Constraint(model.T, rule=TempconstrDHWMin)

def TempconstrDHWMax(model, t):
    return model.T_TES[t] <= model.T_TES_UB

model.TempconstrDHWMax = pyo.Constraint(model.T, rule=TempconstrDHWMax)

def lossDHW(model, t):
    return model.epsilon_TES[t] == model.K_TES * (model.T_TES[t] - model.T_out[t])

model.lossDHW = pyo.Constraint(model.T, rule=lossDHW)

def TempEqvDHW(model):
    return model.T_TES[model.t_end] == model.T_TES_init

model.TempEqvDHW = pyo.Constraint(rule=TempEqvDHW)
      
# Thermal energy storage constraints:
def energyConstrTES(model, t):
    if t == 1:
        return model.Q_TES[t] == model.Q_TES_init + model.q_DHW[t] * model.delta_t - model.q_TES_d[t] * model.delta_t - model.epsilon_TES[t] * model.delta_t
    else:
        return model.Q_TES[t] == model.Q_TES[t-1] + model.q_DHW[t] * model.delta_t - model.q_TES_d[t] * model.delta_t - model.epsilon_TES[t] * model.delta_t
    
model.energyConstrTES = pyo.Constraint(model.T, rule=energyConstrTES)

def energyTESBalance(model):
    return model.Q_TES[model.t_end] == model.Q_TES_init

model.energyTESBalance = pyo.Constraint(rule=energyTESBalance)

def energyMinTES(model, t):
    return model.Q_TES_min <= model.Q_TES[t]

model.energyMinTES = pyo.Constraint(model.T, rule=energyMinTES)

def energyMaxTES(model, t):
    return model.Q_TES[t] <= model.Q_TES_max

model.energyMaxTES = pyo.Constraint(model.T, rule=energyMaxTES)

def powerDischargeMinTES(model, t):
    return model.q_TES_d_min <= model.q_TES_d[t]

model.powerDischargeMinTES = pyo.Constraint(model.T, rule=powerDischargeMinTES)

def powerDischargeMaxTES(model, t):
    return model.q_TES_d[t] <= model.q_TES_d_max

model.powerDischargeMaxTES = pyo.Constraint(model.T, rule=powerDischargeMaxTES)

def TempTES(model, t):
    return model.Q_TES[t] == (model.V_TES * model.rho_TES * model.c_TES * (model.T_TES[t] - model.T_init[t]))/3.6e6

model.TempTES = pyo.Constraint(model.T, rule=TempTES)

# Power balance constraints:
def powerBalance(model, t):
    return model.p_HP[t] * model.delta_t + model.d_ele[t] + model.p_export[t] == model.p_pv[t] * model.delta_t + model.p_import[t]

model.powerBalance = pyo.Constraint(model.T, rule=powerBalance)

# Heat balance constraints:
def heatBalance(model, t):
    return model.d_heat[t] == (model.q_SH[t] + model.q_TES_d[t]) * model.delta_t

model.heatBalance = pyo.Constraint(model.T, rule=heatBalance)

############################################################################
######################## Setting objective function: ####################### 
############################################################################

def ObjectiveFuction(model):
    total = 0 
    for t in model.T:
        total += model.pi_import[t] * model.p_import[t] * model.delta_t - model.pi_export[t] * model.p_export[t] * model.delta_t
    return total

model.obj = pyo.Objective(rule=ObjectiveFuction, sense=pyo.minimize)