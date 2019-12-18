process poly = [ 
    Poisson = 0.3		%Poisson's Ratio = 0.3
    thermcond = 2.33		%Thermal conductivity Si = 2.33e-6/C
    viscosity = 1.78e-5		%Viscosity (of air) = 1,78e-5
    fluid = 2e-6		% Between the device and the substrate.
    density = 2300		%Material density = 2300 kg/m^3
    Youngsmodulus = 165e9	%Young's modulus = 1.65e11 N/m^2
    permittivity = 8.854e-12	%permittivity: C^2/(uN.um^2)=(C.s)^2/kg.um^3;
    sheetresistance = 20	%Poly-Si sheet resistance [ohm/square]
    stress = 0
    straingradient=0
    thermalexpansion=0
    ambienttemperature=0
    
    thermal_coefficient_of_resistance = 2500e-6 %units [1/Kelvin]
    thermal_coefficient_of_expansion = 2.33e-6 %units [1/Kelvin]
    sheet_resistance = 1500 %units [Seiens/meter]. resistance = length/(sheet_resistance * cross_sectional_area).
    thermal_conductivity = 148 %units [Watt/Kelvin/meter]
    electrical_conductivity = 1500 %units [Seimens/meter]
    Youngs_modulus = 160e9 %units [Pascals] [force/area]
    
    ox = 0
    oy = 0
    oz = 0
]

process p1 : poly = [ 
    h = 2e-6			%Layer height of mcnc poly1 = 2e-6 meters
]

process p2 : poly = [
    h = 1.5e-6			%Layer height of mcnc poly2 = 1.5e-6 meters
]

process d2 : poly = [
    fluid = 0.75e-6		%Fluid layer thickness p1:p2 = 0.75e-6 meters.
]

process cmos : poly = [ 
    linewidth = 2
    nodewidth = 6
    ox1=0 
    oy1=0 
    oz1=0 
    ox2=0 
    oy2=0 
    oz2=0 
    L1=0 
    L2=0 
    ox=0
    oy=0
    oz=0
]