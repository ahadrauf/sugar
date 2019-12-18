function trapezoidalbeam

E_ = 170e9;
h = 2e-6;
D = 2300;
len = 100e-6;
L = 200e-6;
V = h*len*len;
m = D*V;

i=0;
if 1 
for wb=0:0.01e-6:0.1e-6;
   wa=2e-6;
   i = i + 1;
   Wa(i) = wa;
   Wb(i) = wb;
   Domega_wb(i) = 1/8.*E_.*h.*(3.*wb.^2+8.*wb.*wa+6.*wa.^2)/(sqrt(E_.*h.*(wb.^3+4.*wb.^2.*wa+6.*wb.*wa.^2+4.*wa.^3)/(L.^3.*m)).*L.^3.*m);
   omega_wb(i) = 1/4.*sqrt(E_.*h.*(wb.^3+4.*wb.^2.*wa+6.*wb.*wa.^2+4.*wa.^3)/(L.^3.*m));
end
figure(1);
plot(Wb,Domega_wb); grid on;
figure(2);
plot(Wb,omega_wb); grid on;
end
   
if 0
for wa=0:0.01e-6:10e-6;
   wb=0e-6;
   i = i + 1;
   Wa(i) = wa;
   Wb(i) = wb;
   Domega_wa(i) = 1/8.*E_.*h.*(4.*wb.^2+12.*wb.*wa+12.*wa.^2)/(sqrt(E_.*h.*(wb.^3+4.*wb.^2.*wa+6.*wb.*wa.^2+4.*wa.^3)/(L.^3.*m)).*L.^3.*m);
   omega_wa(i) = 1/4.*sqrt(E_.*h.*(wb.^3+4.*wb.^2.*wa+6.*wb.*wa.^2+4.*wa.^3)/(L.^3.*m));
end
figure(1);
plot(Wa,Domega_wa); grid on;
figure(2);
plot(Wa,omega_wa); grid on;
end



if 0
j=0;
for wb=0e-6:0.01e-6:0.1e-6;
   i = i + 1;
   for wa=05e-6:0.1e-6:10e-6;
      j = j + 1;      
      Wb(i) = wb;
      Wa(j) = wa;
      Domega_wa(i,j) = 1/8.*E_.*h.*(4.*wb.^2+12.*wb.*wa+12.*wa.^2)/(sqrt(E_.*h.*(wb.^3+4.*wb.^2.*wa+6.*wb.*wa.^2+4.*wa.^3)/(L.^3.*m)).*L.^3.*m);
      omega(i,j) = 1/4.*sqrt(E_.*h.*(wb.^3+4.*wb.^2.*wa+6.*wb.*wa.^2+4.*wa.^3)/(L.^3.*m));
   end
   j=0;
end 

figure(1);
colormap('jet');
shading interp;
surf(Wa, Wb,omega); grid on;
figure(2);
colormap('jet');
shading interp;
surf(Wa,Wb,Domega_wa); grid on; 
end