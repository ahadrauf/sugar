%generalprocess

process fabrication = 
[ 
%    Youngsmodulus = 165e9	%Poly-Si, Young's modulus N/m^2

%    Youngsmodulus = 129.5e9  %[100] SCS, Young's modulus N/m^2
%    shearmodulus  =  79.0e9  %[100] SCS, shear modulus N/m^2

    Youngsmodulus = (129.5e9+186.5e9)/2  %[100]&[110] SCS, Young's modulus N/m^2
    shearmodulus  =  (79.0e9+57.5e9)/2  %[100]&[110] SCS, shear modulus N/m^2
    
%    Youngsmodulus = 168.0e9	%[110] SCS, Young's modulus N/m^2
%    shearmodulus  =  61.7e9  %[110] SCS, shear modulus N/m^2
    
%    Youngsmodulus = 186.5e9	%[111] SCS, Young's modulus N/m^2
%    shearmodulus  =  57.5e9  %[111] SCS, shear modulus N/m^2
    
    Poisson = 0.3		%Poisson's Ratio = 0.3
    thermcond = 2.33		%Thermal conductivity Si = 2.33e-6/C
    viscosity = 1.78e-5		%Viscosity (of air) = 1,78e-5
    fluid = 2e-6		% Between the device and the substrate.
    density = 2300		%Material density = 2300 kg/m^3
    permittivity = 8.854e-12	%permittivity: C^2/(uN.um^2)=(C.s)^2/kg.um^3;
    sheetresistance = 20	%Poly-Si sheet resistance [ohm/square]
    stress = 0
    straingradient=0
    thermalexpansion=0
    ambienttemperature=0
]

process layer1 : fabrication = 
[ 
    h = 2e-6			%Layer height 
    layer = 'p1'
]

process layer2 : fabrication = 
[
    h = 1.5e-6			%Layer height 
    layer = 'p2'
]

process d2 : fabrication = 
[
    fluid = 0.75e-6		%Fluid layer thickness 
]


