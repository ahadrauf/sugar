W=18.858;
E=0.02;
WW=355.31;
F=4.2336e-3 * 600;
C2 = -1/2*W*F*(2*E^2*W^2-2*E*W*sqrt(E^2*W^2-WW)+W^2-WW)/(sqrt(E^2*W^2-WW)*(4*E^2*W^4+WW^2-2*WW*W^2+W^4));
C1 = 1/2*(2*E^2*W^2+2*E*W*sqrt(E^2*W^2-WW)-WW+W^2)*W*F/(sqrt(E^2*W^2-WW)*(4*E^2*W^4+WW^2-2*WW*W^2+W^4));


t=0;
i=0;
figure(1);
clf;

for t = 0:0.01:200
   i=i+1;
   Q(i)=1/2*F*(4*cos(W*t)*W^2*E*sqrt(E^2*W^2-WW)-2*sin(W*t)*E^2*W^2*sqrt(E^2*W^2-WW)+2*sin(W*t)*(E^2*W^2-WW)^(3/2)+2*sin(W*t)*sqrt(E^2*W^2-WW)*W^2)/(sqrt(E^2*W^2-WW)*(-2*E^2*W^2+2*E*W*sqrt(E^2*W^2-WW)+WW-W^2)*(2*E^2*W^2+2*E*W*sqrt(E^2*W^2-WW)-WW+W^2))+C1*exp(-(E*W-sqrt(E^2*W^2-WW))*t)+C2*exp(-(E*W+sqrt(E^2*W^2-WW))*t);
   Qdot=1/2*F*(-4*sin(W*t)*W^3*E*sqrt(E^2*W^2-WW)-2*cos(W*t)*W^3*E^2*sqrt(E^2*W^2-WW)+2*cos(W*t)*W*(E^2*W^2-WW)^(3/2)+2*cos(W*t)*W^3*sqrt(E^2*W^2-WW))/(sqrt(E^2*W^2-WW)*(-2*E^2*W^2+2*E*W*sqrt(E^2*W^2-WW)+WW-W^2)*(2*E^2*W^2+2*E*W*sqrt(E^2*W^2-WW)-WW+W^2))+C1*(-E*W+sqrt(E^2*W^2-WW))*exp(-(E*W-sqrt(E^2*W^2-WW))*t)+C2*(-E*W-sqrt(E^2*W^2-WW))*exp(-(E*W+sqrt(E^2*W^2-WW))*t);
   T(i)=t;
end
plot(T,Q);
grid on;