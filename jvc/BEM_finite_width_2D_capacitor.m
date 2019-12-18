% [sigma,S,V]=BEM_finite_width_2D_capacitor(2e-6,50,10,100e-6);
function [sigma,S,V,Q_bea,C_bea]=BEM_finite_width_2D_capacitor(gap,voltage,number_of_elements_per_plate,width,thickness)
N2=2*number_of_elements_per_plate;
N1=number_of_elements_per_plate;
element_width=width/N1;
a=element_width/2;
S=zeros(N2,N2);
V=zeros(N2,1);x=V;y=V;X=V;Y=V;

element_width*N1

for i=1:N1
    x(i)=-width/2-element_width/2+element_width*i;
    y(i)=gap/2;
    V(i)=voltage/2;
end
for i=1:N1
    x(N1+i)=-width/2-element_width/2+element_width*i;
    y(N1+i)=-gap/2;
    V(N1+i)=-voltage/2;
end
X=x;Y=y;
for j=1:N2
    for i=1:N2
        sx=x(j)-X(i);
        sy=y(j)-Y(i);
        S(j,i)=BEM_S_xzplane(sx,sy,a);
    end
end

e0=8.854187817e-12; %permittivity of free space.
sigma=S\V;
area = width*thickness;
Q_bea=sum(sigma(1:N1)) * area;
C_bea=Q_bea/voltage;
C_bea
C_inf=e0*area/gap
Q_bea
Q_inf=C_inf*voltage


;i=1:N1;figure(1);plot(i*width,sigma(i),'-*b')
grid on; xlabel('length along 2D plate');ylabel('surface charge density');title('Boundary Element Method');
hold on; 
plot(i*width,Q_inf*i./i/area,'-r');

dq_bea=sigma(N1/2)*thickness*element_width
dq_ing= ((e0*thickness*element_width)/gap) * voltage
