taylor expand about x   Cx+s2*x / (  (DC1+s2*x) / (glo + (-glo * (n + (DC2+s2*x) / (DC1+s2*x))/(1 + (DC2+s2*x)/(DC))  )

taylor expand about dC   Cx+s2*dC / (  (DC1+s2*dC) / (glo + (-glo * (n + (DC2+s2*dC) / (DC1+s2*dC))/(1 + (DC2+s2*dC)/(DC))  )

psi = DC1 / (glo + Dg)

Dg = -glo * (n + DC2/DC1)/(1 + DC2/DC1)

x = Cx / psi

%%%%%%%%%%%%%%%

x = Cx / (DC1 / (glo + (-glo * (n + DC2/DC1)/(1 + DC2/DC1))))  
x = Cx / ((CD1 + s2*dC) / (glo + (-glo * (n + (CD2 + s2*dC)/(CD1 + s2*dC))/(1 + (CD2 + s2*dC)/(CD1 + s2*dC)))))  

xpdx = Cx / ((CD1 + s2*x) / (glo + (-glo * (n + (CD2 + s2*dC)/(CD1 + s2*x))/(1 + (CD2 + s2*x)/(CD1 + s2*x)))))  


DC1 = (C1 + s2*x)
DC2 = (C2 + s2*x)

taylor series  Cx / ((CD1 + s2*x) / (glo + (-glo * (n + (CD2 + s2*dC)/(CD1 + s2*x))/(1 + (CD2 + s2*x)/(CD1 + s2*x)))))
taylor series  C / ((C1 + s2*x) / (g + (-g * (n + (C2 + s2*x)/(C1 + s2*x))/(1 + (C2 + s2*x)/(C1 + s2*x)))))

%%%%%%%%%%%%%%%%%

Dg = (-glo * (n + (C2 + s2*x)/(C1 + s2*x))/(1 + (C2 + s2*x)/(C1 + s2*x)))
F = 1/2 * ( (C1 + s2*dC) / (glo + (-glo * (n + (C2 + s2*x)/(C1 + s2*x))/(1 + (C2 + s2*x)/(C1 + s2*x))))) * (V+s2*y)^2

%%%%%%%%%%%%%

s2=sqrt(2); N=46*8*2; e0=8.854e-12; h=20e-6; g=2e-6; x1=10e-6; x2=16e-6; C1=N*e0*h*x1/g; C2=N*e0*h*x2/g; Cx=N*e0*h*5e-6/g; n=1.6; dx_dC=2*s2*g*(n-1)*Cx/(C1+C2)^2, V=10; dF_dC = s2*V^2/(g*(n-1)), dF_dV=s2*(C1+C2)/(g*(n-1))

