import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from rocketcea.cea_obj import CEA_Obj

from propellent_specs_from_pressure import PropellentSpecs

# Conversions
psi_to_pa = 6894.757293168
inches_to_meters = 0.0254

# Constants
tank_id_inches = 7.75 #in
tank_height_inches = 10.5 * 12 # ft --> in
tank_volume_inch_cube = np.pi*np.square((tank_id_inches/2))*tank_height_inches
tank_volume_meter_cube =  tank_volume_inch_cube * (inches_to_meters ** 3)

OF_Ratio = 5.0
prop_tank_to_cc_pressure = 0.6


def pressure_drop_over_sf_burn(propellent_tank_start_pressure, graph=False):
    # Extracting Important LSF5 Data
    # LSF5 Data CSV
    LSF5_df = pd.read_csv("LSF5_data.csv")
    time = LSF5_df["time"]
    prop_tank_pressure_data_psi = LSF5_df["[OPT-102] Ox Tank"]
    cc_pressure_data_psi = LSF5_df["[CCPT-101] Combustion Chamber PT1 (psi)"]

    sf_tank_start_pressure_psi = prop_tank_pressure_data_psi[0]
    
    # Scaling SF Pressure Curve to match simulated starting pressure 
    pressure_scaling_factor_psi = propellent_tank_start_pressure / sf_tank_start_pressure_psi

    prop_tank_pressure_data_psi *= pressure_scaling_factor_psi
    burn_tank_start_pressure_psi = sf_tank_start_pressure_psi * pressure_scaling_factor_psi
    cc_pressure_data_psi *= pressure_scaling_factor_psi
    
    theoretical_cc_pressure_data_psi = prop_tank_pressure_data_psi * prop_tank_to_cc_pressure  # Theoretical CC Pressure is 0.6 Tank pressure
    
    # Finding pressure at burn end time (time when cc pressure drops substatnailly)
    burn_end_time = 2930.0
    diffs = np.abs(time.values - burn_end_time)
    burn_end_index = int(np.argmin(diffs))
    # print("Nearest index:", burn_end_index, "time value:", time.iloc[burn_end_index])

    burn_tank_end_pressure_psi = prop_tank_pressure_data_psi[burn_end_index]

    print(f"Burn Start Pressure: {burn_tank_start_pressure_psi} psi")
    print(f"Burn end Pressure: {burn_tank_end_pressure_psi} psi")

    # Turnicate theoretical cc data when propellents run out and make a copt to graph (Dosen't lag prop.tank during burn like actual CC data) -- Check that assumption is correct May be source of error for almost constant isp???
    graph_theoretical_cc_pressure_data_psi = theoretical_cc_pressure_data_psi
    graph_theoretical_cc_pressure_data_psi[burn_end_index:] = np.ma.masked
    theoretical_cc_pressure_data_psi = theoretical_cc_pressure_data_psi[:burn_end_index]

    # Similarly tunicate tank pressure data for weighted avg
    graph_prop_tank_pressure_data_psi = prop_tank_pressure_data_psi
    graph_prop_tank_pressure_data_psi[burn_end_index:] = np.ma.masked
    prop_tank_pressure_data_psi = prop_tank_pressure_data_psi[:burn_end_index]
    
    # Propellent Tank Pressure Drop Over Burn
    burn_propellent_tank_pressure_drop = (burn_tank_start_pressure_psi - burn_tank_end_pressure_psi ) / burn_tank_start_pressure_psi
    propellent_tank_residual_pressure =  1 - burn_propellent_tank_pressure_drop

    print(f"Propellent tank pressure drop over burn is {burn_propellent_tank_pressure_drop*100} % ")
    # print(f"Sanity Check -- Calculated end pressure from Pressure drop: Prop Tank End pressure = {burn_start_pressure_psi*propellent_tank_residual_pressure} psi")

    # Graphing Prop Tank and CC Pressure 
    if graph:
        plt.plot(time, graph_prop_tank_pressure_data_psi, label="Prop. Tank")
        plt.plot(time, cc_pressure_data_psi, label="Actual CC")
        plt.plot(time, graph_theoretical_cc_pressure_data_psi, label="Theoretical CC")
        plt.xlabel("time")
        plt.ylabel("Pressure (psi)")
        plt.legend()
        plt.show()

    return prop_tank_pressure_data_psi, theoretical_cc_pressure_data_psi, burn_tank_end_pressure_psi 


def isp_calc(cc_pressure_data_psi, OF_Ratio):
    # Calculates Avg. ISP Over Burn -- Assumes CC Pressure follows tank pressure exactlyu *0.6 (No pressure build up in first secs)
    cea = CEA_Obj(oxName="N2O", fuelName='Ethanol')
    
    cc_start_pressure_psi = cc_pressure_data_psi[0]
    ambient_pressure_psi = 14.7 # sea level, or 15.04 at Timmins, ON
    
    # Find expansion ration using call to CEA -- Only dependo on CC Pressure at Start of burn (Full Expansion at Sea Level)
    expansion_ratio = cea.get_eps_at_PcOvPe(Pc=cc_start_pressure_psi, MR=OF_Ratio, PcOvPe=cc_start_pressure_psi/ambient_pressure_psi)
    print(f"esp={expansion_ratio}")

    isp_data = []
    for cc_pressure_psi in theoretical_cc_pressure_data_psi:

        # Get rid of all commented below once sure eps from above is accurate
        # gamma = ... #Whatever cea call for gamma

        # # Exit Mach -- Rearagned From Mach & isentropic pressure relation (since stagnation pressure ~= cc pressure)
        # Me = np.sqrt((2.0/(gamma - 1)) * (((cc_pressure_psi/ambient_pressure_psi)**((gamma - 1)/gamma)) - 1))

        # # Expansion Ratio -- from Nozzle Area Ratio for isentropic flow eqn
        # sqrt_coeff_1 = 2 / (gamma + 1)
        # sqrt_coeff_2 = 1 + ((gamma - 1) / 2)* Me**2
        # sqrt_power = (gamma + 1) / (gamma - 1)
        # expansion_ratio = (1/Me) * np.sqrt((sqrt_coeff_1*sqrt_coeff_2) ** sqrt_power)
        # print("EPS=",expansion_ratio)

        # CEA call to get instant isp at current pressure
        instant_isp = cea.isp = cea.get_Isp(Pc=cc_pressure_psi,
                        MR=OF_Ratio,
                        eps=expansion_ratio)
        
        isp_data.append(instant_isp)
        # print(f"CC Pressure={cc_pressure_psi} psi, instant isp={instant_isp} s") --- Note isp dosen't change much over burn even when pressure drops over 150 psi. SOmthing may be wrong -- Might need to change eps???
    
    isp_avg = np.average(isp_data)
    isp_weighted_avg = np.average(isp_data, weights=prop_tank_pressure_data_psi)

    print(f"Avg isp: {isp_avg} s")
    print(f"Weighted Avg isp: {isp_weighted_avg} s")
        
    return isp


def residual_vapour_mass(burn_tank_end_pressure_psi, prop_tank_volume_meter_cube):
    # Could get Ullage at every time step to try to get an m dot from piston being pushed more and more --> Help to get a weighted avg of ISP
    pressure_pa = burn_tank_end_pressure_psi * psi_to_pa

    ox_vapour_density = CP.PropsSI('D', 'P', pressure_pa, 'Q', 1, 'N2O')
    ox_vapour_mass = ox_vapour_density * prop_tank_volume_meter_cube 

    return ox_vapour_mass

prop_tank_start_pressure_psi = 800

prop_tank_pressure_data_psi, theoretical_cc_pressure_data_psi, burn_tank_end_pressure_psi = pressure_drop_over_sf_burn(prop_tank_start_pressure_psi, graph=True)
isp = isp_calc(theoretical_cc_pressure_data_psi, OF_Ratio)

reisidual_propellent_mass = residual_vapour_mass(burn_tank_end_pressure_psi, tank_volume_meter_cube)
print(f"\nResidual Propellent (N2O Vapour) Mass: {reisidual_propellent_mass} kg")


def calc_delta_v(avg_isp, dry_mass, reisidual_propellent_mass):
    pass

# Next Steps:
# Find correct call to get gamma
# get avg isp and weighted avg isp wrt pressure over burn
# Finally, get dv using dry_mass = 135 lb + residual prop mass, wet_mass = 135 + inital propellent (total) mass (incl. ullage) and avg. isp

# Then, graph dv v pressure from 400 psi to 800 psi to find max dv and Graph Results.


 





