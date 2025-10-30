import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from rocketcea.cea_obj import CEA_Obj

from propellent_specs_from_pressure import PropellentSpecs

# Conversions
psi_to_pa = 6894.757293168
inches_to_meters = 0.0254
prop_tank_to_cc_pressure = 0.6

def isp_calc(cc_pressure_psi, OF_Ratio):
    cea = CEA_Obj(oxName="N2O", fuelName='Ethanol')

    ambient_pressure_psi = 14.7 # sea level, or 15.04 at Timmins, ON

    # # Find Expansion Ratio using Gamma -- Needed for ISP
    gamma = cea.get_gamma(Pc=cc_pressure_psi, MR=OF_Ratio)

    # # Exit Mach -- Rearagned From Mach & isentropic pressure relation (since stagnation pressure ~= cc pressure)
    Me = np.sqrt((2.0/(gamma - 1)) * (((cc_pressure_psi/ambient_pressure_psi)**((gamma - 1)/gamma)) - 1))

    # Expansion Ratio -- from Nozzle Area Ratio for isentropic flow eqn
    sqrt_coeff_1 = 2 / (gamma + 1)
    sqrt_coeff_2 = 1 + ((gamma - 1) / 2)* Me**2
    sqrt_power = (gamma + 1) / (gamma - 1)
    expansion_ratio = (1/Me) * np.sqrt((sqrt_coeff_1*sqrt_coeff_2) ** sqrt_power)
    print("EPS=",expansion_ratio)

    isp = cea.get_Isp(Pc=cc_pressure_psi, MR=OF_Ratio, eps=expansion_ratio)

    return isp

isp_calc(460*0.6, 5)

def calc_delta_v(avg_isp, dry_mass, reisidual_propellent_mass):
    pass

# Next Steps:
# Find correct call to get gamma
# get avg isp and weighted avg isp wrt pressure over burn
# Finally, get dv using dry_mass = 135 lb + residual prop mass, wet_mass = 135 + inital propellent (total) mass (incl. ullage) and avg. isp

# Then, graph dv v pressure from 400 psi to 800 psi to find max dv and Graph Results.


 





