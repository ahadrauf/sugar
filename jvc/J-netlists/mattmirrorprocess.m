process poly=[
   Poisson=0.3 			 %Poisson's Ratio = 0.3
   thermcond=2.33 		 %Thermal conductivity Si = 2.33e-6/C
   viscosity=1.78e-5 	 %Viscosity (of air) = 1,78e-5
   fluid=2e-6 				 %Between the device and the substrate.
   density=2300 			 %Material density = 2300 kg/m^3
   Youngsmodulus=185e9	 %Young's modulus = 1.65e11 N/m^2
   permittivity=8.854e-12%permittivity: C^2/(uN.um^2)=(C.s)^2/kg.um^3;
   sheetresistance=20 	 %Poly-Si sheet resistance [ohm/square]
   stress=12e6				 %Residual planar stress
   straingradient=0		 %Residual strain gradient
   thermalexpansion=0 	 %Coefficient of thermal expansion
   ambienttemperature=0	 %Ambient Temperature
]
 
process p1:poly=[ 
   h=40e-6 %Layer height of mcnc poly1 = 2e-6 meters
]

process p2:poly=[
   h=4e-6 %Layer height of mcnc poly2 = 1.5e-6 meters
]

process d2:poly=[
   fluid=0.75e-6 %Fluid layer thickness p1:p2 = 0.75e-6 meters
]

