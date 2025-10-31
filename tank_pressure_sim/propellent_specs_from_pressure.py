"""
Dynamically find propellents mass and each propellents tank volume from tank pressure

Design/Independant Variables: Pressure, ullage_by_volume (fraction by volume)

Assumptions: Constant Ullage by Volume
"""

import CoolProp.CoolProp as CP
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Conversions
psi_to_pa = 6894.757293168
inches_to_meters = 0.0254
prop_tank_to_cc_pressure = 0.6


def oxTankVolume(ullage_by_volume, OF_Ratio, tank_volume_meter_cube, ox_liquid_density, fuel_density):
    
    ox_volume = (OF_Ratio * fuel_density * tank_volume_meter_cube) / (((1 - ullage_by_volume) * ox_liquid_density) + OF_Ratio * fuel_density)

    return ox_volume


def PropellentSpecs(pressure_psi, ullage_by_volume, OF_Ratio, tank_volume_meter_cube):
    # Oxizizer Mass
    # --> Mass = volume x density
    # --> Density: updates with pressure from changing tank temp.
    pressure_pa = pressure_psi * psi_to_pa

    ox_liquid_density = CP.PropsSI('D', 'P', pressure_pa, 'Q', 0, 'N2O')
    ox_vapour_density = CP.PropsSI('D', 'P', pressure_pa, 'Q', 1, 'N2O')
    fuel_density = CP.PropsSI('D', 'P', pressure_pa, 'Q', 0, 'C2H6O')

    # print(f'\nliquid ox density: {ox_liquid_density} kg/m^3')
    # print(f'vapour ox density: {ox_vapour_density} kg/m^3')
    # print(f'fuel density: {fuel_density} kg/m^3')
    
    ox_tank_volume_meter_cube = oxTankVolume(ullage_by_volume, OF_Ratio, tank_volume_meter_cube, ox_liquid_density, fuel_density)
    fuel_tank_volume_meter_cube = tank_volume_meter_cube - ox_tank_volume_meter_cube
    
    ox_mass = ox_liquid_density * ((1 - ullage_by_volume)* ox_tank_volume_meter_cube)
    fuel_mass = fuel_density * fuel_tank_volume_meter_cube
    ox_vapour_mass = (ullage_by_volume * ox_tank_volume_meter_cube) * ox_vapour_density

    return ox_mass, fuel_mass, ox_vapour_mass, ox_tank_volume_meter_cube, fuel_tank_volume_meter_cube


def propellent_mass_by_pressure_graph(ullage_by_volume, OF_Ratio, tank_volume_meter_cube):
    # Graph -- Propellent Masses vs Pressure 
    pressure_min_psi, pressure_max_psi, pressure_step_psi = 400, 800, 50

    psi_pressures = np.arange(pressure_min_psi, pressure_max_psi, pressure_step_psi) 

    fuel_masses = []
    ox_liquid_masses = []
    ox_vapour_masses = []
    
    ox_tank_volumes_inches_cube = []
    fuel_tank_volumes_inches_cube = []

    for P_psi in psi_pressures:
        ox_mass, fuel_mass, ox_vapour_mass, ox_tank_volume_meter_cube, fuel_tank_volume_meter_cube = PropellentSpecs(P_psi, ullage_by_volume, OF_Ratio, tank_volume_meter_cube)
       
        fuel_masses.append(fuel_mass)
        ox_liquid_masses.append(ox_mass)
        ox_vapour_masses.append(ox_vapour_mass)
        
        ox_tank_volumes_inches_cube.append(ox_tank_volume_meter_cube / (inches_to_meters ** 3))
        fuel_tank_volumes_inches_cube.append(fuel_tank_volume_meter_cube / (inches_to_meters ** 3))

    fuel_masses = np.array(fuel_masses)
    ox_liquid_masses = np.array(ox_liquid_masses)
    ox_vapour_masses = np.array(ox_vapour_masses)
    ox_tank_volumes_inches_cube = np.array(ox_tank_volumes_inches_cube)
    fuel_tank_volumes_inches_cube = np.array(fuel_tank_volumes_inches_cube)

    print('From Graphs:')
    print(f'Fuel Mass -- Max: {max(fuel_masses)} kg, Min: {min(fuel_masses)} kg')
    print(f'Liquid Ox Mass -- Max: {max(ox_liquid_masses)} kg, Min: {min(ox_liquid_masses)} kg')
    print(f'Vapour Ox Mass -- Max: {max(ox_vapour_masses)} kg, Min: {min(ox_vapour_masses)} kg \n')

    print(f'Ox Tank Volume -- Max: {max(ox_tank_volumes_inches_cube)} in^3, Min: {min(ox_tank_volumes_inches_cube)} in^3')
    print(f'Fuel Tank Volume -- Max: {max(fuel_tank_volumes_inches_cube)} in^3, Min: {min(fuel_tank_volumes_inches_cube)} in^3')

    plt.figure(1)
    plt.plot(psi_pressures, ox_liquid_masses + fuel_masses,   label='Propellent mass (kg)')
    plt.plot(psi_pressures, fuel_masses, label='Fuel mass (kg)')
    plt.plot(psi_pressures, ox_liquid_masses,   label='Liquid Ox mass (kg)')
    plt.plot(psi_pressures, ox_vapour_masses,   label='Vapour Ox mass (kg)')
    plt.xlabel('Tank Pressure (psi)')
    plt.ylabel('Mass (kg)')
    plt.legend(); plt.grid(True)

    plt.figure(2)
    plt.plot(psi_pressures, ox_tank_volumes_inches_cube,   label='Ox Tank Volume (in^3)')
    plt.plot(psi_pressures, fuel_tank_volumes_inches_cube, label='Fuel Tank Volume (in^3)')
    plt.xlabel('Tank Pressure (psi)')
    plt.ylabel('Propellent Tank Volumes (in^3)')
    plt.legend(); plt.grid(True); 
    
    plt.show()

def main():
    
    # Constants
    tank_id_inches = 7.75 #in
    tank_height_inches = 10.5 * 12 # ft --> in
    tank_volume_inch_cube = np.pi*np.square((tank_id_inches/2))*tank_height_inches
    tank_volume_meter_cube =  tank_volume_inch_cube * (inches_to_meters ** 3)
    
    # Design Variables
    pressure_psi = 650 
    ullage_by_volume = 0.12
    OF_Ratio = 5.0
    
    # Propellent Mass and Tank Specs
    ox_mass, fuel_mass, ox_vapour_mass, ox_tank_volume_meter_cube, fuel_tank_volume_meter_cube = PropellentSpecs(pressure_psi, ullage_by_volume, OF_Ratio, tank_volume_meter_cube)

    ox_tank_volume_inch_cube = ox_tank_volume_meter_cube / (inches_to_meters ** 3)
    fuel_tank_volume_inch_cube = fuel_tank_volume_meter_cube / (inches_to_meters ** 3)

    print(f'Test @ {pressure_psi} psi:')
    print(f'Pressure in psi: {pressure_psi}')
    print(f'Pressure in pa: {pressure_psi * psi_to_pa}')

    print(f'\nTank Volume: {tank_volume_inch_cube} in^3')
    print(f'Oxidizer Tank Volume: {ox_tank_volume_inch_cube} in^3')
    print(f'Fuel Tank Volume: {fuel_tank_volume_inch_cube} in^3')

    print(f'\nOxidizer (Liquid) Mass: {ox_mass} kg')
    print(f'Oxidizer (Vapour) Mass: {ox_vapour_mass} kg')
    print(f'Fuel Mass: {fuel_mass} kg')

    print(f'\nSanity Check --  Calculated O/F Ratio: {ox_mass/fuel_mass} \n')


if __name__ == '__main__':
    main()