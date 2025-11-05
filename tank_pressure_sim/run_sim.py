import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from propellent_specs_from_pressure import propellent_specs, propellent_mass_by_pressure_graph
from delta_v_from_pressure import pressure_drop_over_sf_burn, residual_vapour_mass, isp_calc, delta_v_calc

# Conversions & Constants
psi_to_pa = 6894.757293168
inches_to_meters = 0.0254
lbm_to_kg = 0.4536  # Weight is in kg unless otherwise denoted
g_imp = 32.174 #ft/s^2

# Rocket Specs
tank_id_inches = 7.75 #in
tank_height_inches = 10.5 * 12 # ft --> in
tank_volume_inch_cube = np.pi*np.square((tank_id_inches/2))*tank_height_inches
tank_volume_meter_cube =  tank_volume_inch_cube * (inches_to_meters ** 3)
vehicle_mass_lb = 135
vehicle_mass_kg = vehicle_mass_lb * lbm_to_kg
prop_tank_to_cc_pressure = 0.6

# Design Variables
ullage_by_volume = 0.12
prop_tank_to_cc_pressure = 0.6
OF_Ratio = 5.0


# Graph Propellent Masses and Tank Volume from 400 to 800 psi (in incerments of 50 psi) 
def tank_pressure_sim(pressure_min_psi, pressure_max_psi_exclusive, pressure_step_psi, one_graph):

    psi_pressures = np.arange(pressure_min_psi, pressure_max_psi_exclusive, pressure_step_psi) 
    
    # Create empty DataFrame with all specs 
    sim_df = pd.DataFrame(columns=["tank_pressure_psi", "delta_v", "avg_isp", "ox_temp_deg_c", "fuel_mass", "ox_mass", "fuel_tank_volume_inch_cube", "ox_tank_volume_inch_cube"])

    for tank_pressure_psi in psi_pressures:
        # Propellent Mass and Tank Specs
        ox_mass, fuel_mass, ox_vapour_mass, ox_temp_deg_c, ox_tank_volume_meter_cube, fuel_tank_volume_meter_cube = propellent_specs(tank_pressure_psi, ullage_by_volume, OF_Ratio, tank_volume_meter_cube)

        # simulated Propellent and CC Tank Pressure data over burn
        prop_tank_pressure_data_psi, theoretical_cc_pressure_data_psi, burn_tank_end_pressure_psi = pressure_drop_over_sf_burn(tank_pressure_psi, prop_tank_to_cc_pressure, graph=False)
        
        # calculate weighted avg isp over burn
        isp_avg = isp_calc(theoretical_cc_pressure_data_psi, prop_tank_pressure_data_psi, OF_Ratio)

        # residual mass of propellent (oxidizer) after burn
        reisidual_propellent_mass = residual_vapour_mass(burn_tank_end_pressure_psi, tank_volume_meter_cube)

        # masses for delta v calc
        initial_propellent_mass = fuel_mass + ox_mass + ox_vapour_mass
        wet_mass = vehicle_mass_kg + initial_propellent_mass
        dry_mass = vehicle_mass_kg + reisidual_propellent_mass

        dv = delta_v_calc(isp_avg, wet_mass, dry_mass)

        sim_df.loc[len(sim_df)] = [tank_pressure_psi, dv, isp_avg, ox_temp_deg_c, fuel_mass, ox_mass, fuel_tank_volume_meter_cube/(inches_to_meters ** 3), ox_tank_volume_meter_cube/(inches_to_meters ** 3)]
    
    # set pressure to index of df
    sim_df = sim_df.set_index("tank_pressure_psi")

    # Max Values
    optimal_pressure_psi = sim_df["delta_v"].idxmax()
    optimal_data = sim_df.loc[optimal_pressure_psi]
    max_dv = optimal_data["delta_v"]
    isp_optimal_pressure = optimal_data["avg_isp"]
    optimal_temp_deg_c = optimal_data["ox_temp_deg_c"]
    optimal_fuel_mass = optimal_data["fuel_mass"]
    optimal_ox_mass = optimal_data["ox_mass"]
    optimal_fuel_tank_volume_inch_cube = optimal_data["fuel_tank_volume_inch_cube"]
    optimal_ox_tank_volume_inch_cube = optimal_data["ox_tank_volume_inch_cube"]

    # Sim Output
    print("TANK PRESSURE SIMULATION COMPLETE")
    print(f"A maximum delta v of {max_dv} ft/s is achieved at a tank pressure of {optimal_pressure_psi} psi")
    print("Specs at this Pressure:")
    print(f"Isp: {isp_optimal_pressure}")
    print(f"Oxidizer Temp: {optimal_temp_deg_c} °C")
    print(f"Fuel Mass: {optimal_fuel_mass} kg, Liquid Oxidizer Mass: {optimal_ox_mass} kg")
    print(f"Fuel Tank Volume: {optimal_fuel_tank_volume_inch_cube} in^3, Oxidizer Tank Volume: {optimal_ox_tank_volume_inch_cube} in^3 \n")
    print("Full Sim Output")
    print(sim_df)

    if one_graph:
        # full Graph of Simulation 
        plt.plot(sim_df.index, sim_df["delta_v"], label="delta v (ft/s)")
        plt.plot(sim_df.index, sim_df["avg_isp"], label="isp (s)")
        plt.plot(sim_df.index, sim_df["ox_temp_deg_c"], label="Oxidizer temp. (°C)")
        plt.plot(sim_df.index, sim_df["fuel_mass"], label="fuel mass (kg)")
        plt.plot(sim_df.index, sim_df["ox_mass"], label="oxidizer mass (kg)")
        plt.xlabel("Pressure (Psi)")
        plt.ylabel("Data with Respective Value")
        plt.legend()
        plt.show()
    else:
    # Seperate Graph of Simulation
        plt.figure(1)
        plt.title("Tank Pressure Sim: Delta V")
        plt.plot(sim_df.index, sim_df["delta_v"], 'r')
        plt.xlabel("Pressure (Psi)")
        plt.ylabel("delta v (ft/s)")

        plt.figure(2)
        plt.title("Tank Pressure Sim: Isp")
        plt.plot(sim_df.index, sim_df["avg_isp"], 'c')
        plt.xlabel("Pressure (Psi)")
        plt.ylabel("Isp (s)")

        plt.figure(3)
        plt.title("Tank Pressure Sim: Oxidizer Temperature")
        plt.plot(sim_df.index, sim_df["ox_temp_deg_c"], 'm')
        plt.xlabel("Pressure (Psi)")
        plt.ylabel("Oxidizer Temp. (°C)")

        plt.figure(4)
        plt.title("Tank Pressure Sim: Propellent Masses")
        plt.plot(sim_df.index, sim_df["fuel_mass"], 'g', label="fuel mass (kg)")
        plt.plot(sim_df.index, sim_df["ox_mass"], 'b', label="oxidizer mass (kg)")
        plt.xlabel("Pressure (Psi)")
        plt.ylabel("Propellent Mass (kg)")
        plt.legend()
        plt.show()

    return sim_df


# Comparison of our rocketry blowdown system (Piloted by Oxidizer (N2O) Boil Off) vs Standard Blowndown Archetecture (Piloted by 3rd inert Gas, in this case Nitrogen)
def blowdown_systems_comparison(blowdown_system_masses, tank_sim_df):
    tank_pressure_psi = 350
    
    sim_df = pd.DataFrame(columns=["blowdown_system_mass", "delta_v", "avg_isp"])
    
    ox_mass, fuel_mass, ox_vapour_mass, _,  _, _ = propellent_specs(tank_pressure_psi, ullage_by_volume, OF_Ratio, tank_volume_meter_cube)
    prop_tank_pressure_data_psi, theoretical_cc_pressure_data_psi, burn_tank_end_pressure_psi = pressure_drop_over_sf_burn(tank_pressure_psi, prop_tank_to_cc_pressure=1.5, graph=False)
    isp_avg = isp_calc(theoretical_cc_pressure_data_psi, prop_tank_pressure_data_psi, OF_Ratio)
    reisidual_propellent_mass = residual_vapour_mass(burn_tank_end_pressure_psi, tank_volume_meter_cube)

    for blowdown_system_mass in blowdown_system_masses:
        initial_propellent_mass = fuel_mass + ox_mass + ox_vapour_mass
        wet_mass = vehicle_mass_kg + initial_propellent_mass + blowdown_system_mass
        dry_mass = vehicle_mass_kg + reisidual_propellent_mass + blowdown_system_mass
        dv = delta_v_calc(isp_avg, wet_mass, dry_mass)

        sim_df.loc[len(sim_df)] = [blowdown_system_mass, dv, isp_avg]

    # set system masses to index of df
    sim_df = sim_df.set_index("blowdown_system_mass")

    # Max Values
    optimal_mass_kg = sim_df["delta_v"].idxmax()
    optimal_data = sim_df.loc[optimal_mass_kg]
    max_dv = optimal_data["delta_v"]

    # Sim Output
    print("\nARCHITECTURE COMPARISON SIMULATION COMPLETE")
    print(f"A maximum delta v of {max_dv} ft/s is achieved with a Nitrogen Piloted Blowdown Archetecture when the system mass is {optimal_mass_kg} kg")
    print("Full Sim Output")
    print(sim_df)

    # Graphed Comparison
    plt.figure(1)
    plt.title("Tank Architecture Sim: Delta V")
    plt.plot(tank_sim_df.index, tank_sim_df["delta_v"], label="Our Archetecture (psi)")
    plt.plot(sim_df.index, sim_df["delta_v"], label="Nitrogen Piloted (kg)")
    plt.xlabel("Data with Respective Value")
    plt.ylabel("delta v (ft/s)")
    plt.legend()

    plt.figure(2)
    plt.title("Tank Architecture Sim: Isp")
    plt.plot(tank_sim_df.index, tank_sim_df["avg_isp"], label="Our Archetecture (psi)")
    plt.plot(sim_df.index, sim_df["avg_isp"], label="Nitrogen Piloted (kg)")
    plt.xlabel("Data with Respective Value")
    plt.ylabel("Isp (s)")
    plt.legend()
    plt.show()


tank_sim_df = tank_pressure_sim(pressure_min_psi=100, pressure_max_psi_exclusive=850, pressure_step_psi=50, one_graph=False)
blowdown_systems_masses_kg = [5, 10, 15, 20, 25, 30] #kg
blowdown_systems_comparison(blowdown_systems_masses_kg, tank_sim_df)