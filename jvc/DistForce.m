function [F1,M1,F2,M2]=DistForce(equation,L)

%Integration
P=equation;  %units of N/m
psi1=strrep('.*(1-3.*(x./L).^2+2.*(x./L).^3)','L',num2str(L));
psi2=strrep('.*(x.*(1-x./L).^2)','L',num2str(L));
psi3=strrep('.*(3.*(x./L).^2-2.*(x./L).^3)','L',num2str(L));
psi4=strrep('.*(x.^2./L.*(x./L-1))','L',num2str(L)); 

%The forcess
F1=quad8(inline(strcat(P,psi1)),0,L);
M1=quad8(inline(strcat(P,psi2)),0,L);
F2=quad8(inline(strcat(P,psi3)),0,L);
M2=quad8(inline(strcat(P,psi4)),0,L);      

