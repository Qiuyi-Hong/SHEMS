"""
This file contains the code for modelling SHEMS. 
"""

import pyomo.environ as pyo 

def shems_model(t_end):
    model = pyo.AbstractModel(name="SHEMS")

    ############################################################################
    ######################## Setting parameters: ###############################
    ############################################################################
    model.T = pyo.RangeSet(1, t_end, 1)
    model.delta_t = pyo.Param()

    # Electricity demand of the house in kW
    model.d_ele = pyo.Param(model.T, mutable=True)

    # Space heating demand of the house in kW
    model.d_SH = pyo.Param(model.T, mutable=True)

    # DHW demand of the house in kW
    model.d_DHW = pyo.Param(model.T, mutable=True)

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
    model.rho_in = pyo.Param()
    model.V_in = pyo.Param()
    model.c_in = pyo.Param()
    model.T_in_LB = pyo.Param()
    model.T_in_UB = pyo.Param()
    model.K_SH = pyo.Param()
    model.T_in_init = pyo.Param()
    # For DHW:
    model.c_TES = pyo.Param()
    model.T_TES_LB = pyo.Param()
    model.T_TES_UB = pyo.Param()
    model.K_TES = pyo.Param()
    model.T_TES_init = pyo.Param(mutable=True)

    # Thermal energy storage params:
    model.Q_TES_min = pyo.Param()
    model.Q_TES_max = pyo.Param(mutable=True)
    model.Q_TES_init = pyo.Param(mutable=True)
    model.V_TES = pyo.Param(mutable=True)
    model.rho_TES = pyo.Param()
    model.T_inlet = pyo.Param(mutable=True)
    model.T_TES_max = pyo.Param(mutable=True)

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


    # Electricity power imported from the grid in kW
    model.p_import = pyo.Var(model.T, within=pyo.NonNegativeReals)

    # Electricity power exported to the grid in kW
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
            return model.T_in[t] == model.T_in_init + ((model.q_SH[t] - model.d_SH[t] - model.epsilon_SH[t]) * model.delta_t) * 3.6e6 / (model.rho_in * model.V_in * model.c_in)
        else:
            return model.T_in[t] == model.T_in[t-1] + ((model.q_SH[t] - model.d_SH[t] - model.epsilon_SH[t]) * model.delta_t) * 3.6e6 / (model.rho_in * model.V_in * model.c_in)

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

    # def TempEqvSH(model):
    #     return model.T_in[model.t_end] == model.T_in_init

    # model.TempEqvSH = pyo.Constraint(rule=TempEqvSH)

    # For DHW:
    def TempDHW(model, t):
        if t == 1:
            return model.T_TES[t] == model.T_TES_init + ((model.q_DHW[t] - model.d_DHW[t] - model.epsilon_TES[t]) * model.delta_t) * 3.6e6 / (model.rho_TES * model.V_TES * model.c_TES)
        else:
            return model.T_TES[t] == model.T_TES[t-1] + ((model.q_DHW[t] - model.d_DHW[t] - model.epsilon_TES[t]) * model.delta_t) * 3.6e6 / (model.rho_TES * model.V_TES * model.c_TES)

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

    # def TempEqvDHW(model):
    #     return model.T_TES[model.t_end] == model.T_TES_init

    # model.TempEqvDHW = pyo.Constraint(rule=TempEqvDHW)
        
    # Thermal energy storage constraints:
    def energyConstrTES(model, t):
        if t == 1:
            return model.Q_TES[t] == model.Q_TES_init + (model.q_DHW[t] - model.d_DHW[t] - model.epsilon_TES[t]) * model.delta_t
        else:
            return model.Q_TES[t] == model.Q_TES[t-1] + (model.q_DHW[t] - model.d_DHW[t] - model.epsilon_TES[t]) * model.delta_t
        
    model.energyConstrTES = pyo.Constraint(model.T, rule=energyConstrTES)

    # def energyTESBalance(model):
    #     return model.Q_TES[model.t_end] == model.Q_TES_init

    # model.energyTESBalance = pyo.Constraint(rule=energyTESBalance)

    def energyMinTES(model, t):
        return model.Q_TES_min <= model.Q_TES[t]

    model.energyMinTES = pyo.Constraint(model.T, rule=energyMinTES)

    def energyMaxTES(model, t):
        return model.Q_TES[t] <= model.Q_TES_max

    model.energyMaxTES = pyo.Constraint(model.T, rule=energyMaxTES)

    def TempTES(model, t):
        return model.Q_TES[t] == (model.V_TES * model.rho_TES * model.c_TES * (model.T_TES[t] - model.T_inlet))/3.6e6

    model.TempTES = pyo.Constraint(model.T, rule=TempTES)

    # Power balance constraints:
    def powerImportMin(model, t):
        return model.p_import[t] >= 0

    model.powerImportMin = pyo.Constraint(model.T, rule=powerImportMin)

    def powerExportMin(model, t):
        return model.p_export[t] >= 0

    model.powerExportMin = pyo.Constraint(model.T, rule=powerExportMin)

    def powerExportMax(model, t):
        return model.p_export[t] <= model.p_pv[t]

    model.powerExportMax = pyo.Constraint(model.T, rule=powerExportMax)

    def powerBalance(model, t):
        return model.p_HP[t] + model.d_ele[t] + model.p_export[t] == model.p_pv[t] + model.p_import[t]

    model.powerBalance = pyo.Constraint(model.T, rule=powerBalance)

    ############################################################################
    ######################## Setting objective function: ####################### 
    ############################################################################

    def ObjectiveFuction(model):
        total = 0 
        for t in model.T:
            total += (model.pi_import[t] * model.p_import[t] - model.pi_export[t] * model.p_export[t]) * model.delta_t
        return total

    model.obj = pyo.Objective(rule=ObjectiveFuction, sense=pyo.minimize)
    
    return model