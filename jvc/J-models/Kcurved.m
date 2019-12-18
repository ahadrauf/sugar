function [k]=Kcurved(alpha)
w=10e-6;
h=2e-6;
E=165e9;
A=w*h;
EA=E*A;
%Ix=h^3*w/12;
J=w*h^3*(16/3-3.36*h/w*(1-(h/w)^4/12))/16;
Ix=J;
Iy=h*w^3/12;
%alpha=
R=100e-6;
I=Iy;
EI=E*I;
Nu=0.3;
G=E/(2*(1+Nu)); %Shear modulus. GJ = Torsional stiffness 
GIx=G*Ix;

%_______________________________________________________________________________________________________________

k21(1,1)=-2.*EA*GIx*(cos(alpha)*GIx*alpha^2+cos(alpha)*R^2*EA*alpha^2+2.*R^2*EA*cos(alpha)^2-4.*R^2*EA*cos(alpha)-1.*GIx*alpha*sin(alpha)-1.*R^2*EA*alpha*sin(alpha)+2.*R^2*EA)/(R*(4.*R^4*EA^2*alpha*cos(alpha)-4.*R^4*EA^2*sin(alpha)*cos(alpha)+4.*GIx*R^2*EA*sin(alpha)+2.*GIx*alpha*R^2*EA*cos(alpha)^2+4.*GIx*alpha*R^2*EA*cos(alpha)-4.*cos(alpha)*sin(alpha)*GIx*R^2*EA-5.*R^4*EA^2*alpha+GIx^2*alpha^3-1.*GIx^2*alpha+2.*GIx*alpha^3*R^2*EA-6.*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4.*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3));
k21(1,2)=2.*GIx*EA*sin(alpha)*(2.*R^2*EA*cos(alpha)+GIx*alpha^2+R^2*EA*alpha^2-2.*R^2*EA)/(R*(4.*R^4*EA^2*alpha*cos(alpha)-4.*R^4*EA^2*sin(alpha)*cos(alpha)+4.*GIx*R^2*EA*sin(alpha)+2.*GIx*alpha*R^2*EA*cos(alpha)^2+4.*GIx*alpha*R^2*EA*cos(alpha)-4.*cos(alpha)*sin(alpha)*GIx*R^2*EA-5.*R^4*EA^2*alpha+GIx^2*alpha^3-1.*GIx^2*alpha+2.*GIx*alpha^3*R^2*EA-6.*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4.*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3));
k21(1,3)=0;
k21(1,4)=0;
k21(1,5)=0;
k21(1,6)=-2.*EA*GIx*(-1.*cos(alpha)^2*GIx-4.*R^2*EA*cos(alpha)+alpha^2*cos(alpha)*GIx-2.*GIx*alpha*sin(alpha)-2.*R^2*EA*alpha*sin(alpha)+R^2*EA*cos(alpha)^2+cos(alpha)*R^2*EA*alpha^2+GIx+3.*R^2*EA)/(4.*R^4*EA^2*alpha*cos(alpha)-4.*R^4*EA^2*sin(alpha)*cos(alpha)+4.*GIx*R^2*EA*sin(alpha)+2.*GIx*alpha*R^2*EA*cos(alpha)^2+4.*GIx*alpha*R^2*EA*cos(alpha)-4.*cos(alpha)*sin(alpha)*GIx*R^2*EA-5.*R^4*EA^2*alpha+GIx^2*alpha^3-1.*GIx^2*alpha+2.*GIx*alpha^3*R^2*EA-6.*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4.*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3);

k21(2,1)=-2.*GIx*EA*sin(alpha)*(GIx*alpha^2+R^2*EA*alpha^2+2.*R^2*EA*cos(alpha)-2.*R^2*EA)/(R*(4.*R^4*EA^2*alpha*cos(alpha)-4.*R^4*EA^2*sin(alpha)*cos(alpha)+4.*GIx*R^2*EA*sin(alpha)+2.*GIx*alpha*R^2*EA*cos(alpha)^2+4.*GIx*alpha*R^2*EA*cos(alpha)-4.*cos(alpha)*sin(alpha)*GIx*R^2*EA-5.*R^4*EA^2*alpha+GIx^2*alpha^3-1.*GIx^2*alpha+2.*GIx*alpha^3*R^2*EA-6.*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4.*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3));
k21(2,2)=-2.*EA*GIx*(GIx*alpha*sin(alpha)+R^2*EA*alpha*sin(alpha)-2.*R^2*EA+2.*R^2*EA*cos(alpha)^2+alpha^2*cos(alpha)*GIx+cos(alpha)*R^2*EA*alpha^2)/(R*(4.*R^4*EA^2*alpha*cos(alpha)-4.*R^4*EA^2*sin(alpha)*cos(alpha)+4.*GIx*R^2*EA*sin(alpha)+2.*GIx*alpha*R^2*EA*cos(alpha)^2+4.*GIx*alpha*R^2*EA*cos(alpha)-4.*cos(alpha)*sin(alpha)*GIx*R^2*EA-5.*R^4*EA^2*alpha+GIx^2*alpha^3-1.*GIx^2*alpha+2.*GIx*alpha^3*R^2*EA-6.*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4.*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3));
k21(2,3)=0;
k21(2,4)=0;
k21(2,5)=0;
k21(2,6)=-2.*EA*GIx*(-1.*R^2*EA*sin(alpha)-1.*GIx*cos(alpha)*sin(alpha)+alpha*cos(alpha)*R^2*EA+sin(alpha)*GIx*alpha^2-1.*GIx*alpha+sin(alpha)*GIx+alpha*cos(alpha)*GIx+R^2*EA*cos(alpha)*sin(alpha)+sin(alpha)*R^2*EA*alpha^2-1.*R^2*EA*alpha)/(4.*R^4*EA^2*alpha*cos(alpha)-4.*R^4*EA^2*sin(alpha)*cos(alpha)+4.*GIx*R^2*EA*sin(alpha)+2.*GIx*alpha*R^2*EA*cos(alpha)^2+4.*GIx*alpha*R^2*EA*cos(alpha)-4.*cos(alpha)*sin(alpha)*GIx*R^2*EA-5.*R^4*EA^2*alpha+GIx^2*alpha^3-1.*GIx^2*alpha+2.*GIx*alpha^3*R^2*EA-6.*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4.*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3);

k21(3,1)=0;
k21(3,2)=0;
k21(3,3)=-1.*alpha*EI/(R^3*(alpha^2-2.+2.*cos(alpha)));
k21(3,4)=-1.*(alpha-1.*sin(alpha))*EI/(R^2*(alpha^2-2.+2.*cos(alpha)));
k21(3,5)=(cos(alpha)-1.)*EI/(R^2*(alpha^2-2.+2.*cos(alpha)));
k21(3,6)=0;

k21(4,1)=0;
k21(4,2)=0;
k21(4,3)=-1.*(alpha-1.*sin(alpha))*EI/(R^2*(alpha^2-2.+2.*cos(alpha)));
k21(4,4)=-1.*EI*(alpha^2-2.*alpha*sin(alpha)+alpha^2*cos(alpha)+cos(alpha)^2-2.*cos(alpha)+1.)/(R*alpha*(alpha^2-2.+2.*cos(alpha)));
k21(4,5)=EI*(alpha*cos(alpha)-1.*alpha+cos(alpha)*sin(alpha)+sin(alpha)*alpha^2-1.*sin(alpha))/(R*alpha*(alpha^2-2.+2.*cos(alpha)));
k21(4,6)=0;

k21(5,1)=0;
k21(5,2)=0;
k21(5,3)=-1.*EI*(cos(alpha)-1.)/(R^2*(alpha^2-2.+2.*cos(alpha)));
k21(5,4)=-1.*EI*(alpha^2*sin(alpha)-1.*alpha+cos(alpha)*sin(alpha)-1.*sin(alpha)+cos(alpha)*alpha)/(R*alpha*(alpha^2-2.+2.*cos(alpha)));
k21(5,5)=-1.*EI*(-1.+cos(alpha)^2+cos(alpha)*alpha^2)/(R*alpha*(alpha^2-2.+2.*cos(alpha)));
k21(5,6)=0;

k21(6,1)=-2.*EA*GIx*(R^2*EA*cos(alpha)^2+GIx-1.*cos(alpha)^2*GIx-2.*GIx*alpha*sin(alpha)-4.*R^2*EA*cos(alpha)+cos(alpha)*GIx*alpha^2-2.*R^2*EA*alpha*sin(alpha)+cos(alpha)*R^2*EA*alpha^2+3.*R^2*EA)/(4.*R^4*EA^2*alpha*cos(alpha)-4.*R^4*EA^2*sin(alpha)*cos(alpha)+4.*GIx*R^2*EA*sin(alpha)+2.*GIx*alpha*R^2*EA*cos(alpha)^2+4.*GIx*alpha*R^2*EA*cos(alpha)-4.*cos(alpha)*sin(alpha)*GIx*R^2*EA-5.*R^4*EA^2*alpha+GIx^2*alpha^3-1.*GIx^2*alpha+2.*GIx*alpha^3*R^2*EA-6.*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4.*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3);
k21(6,2)=2.*GIx*EA*(sin(alpha)*GIx+sin(alpha)*GIx*alpha^2-1.*GIx*alpha+sin(alpha)*R^2*EA*alpha^2+alpha*cos(alpha)*R^2*EA+R^2*EA*cos(alpha)*sin(alpha)-1.*R^2*EA*sin(alpha)-1.*R^2*EA*alpha+alpha*cos(alpha)*GIx-1.*GIx*cos(alpha)*sin(alpha))/(4.*R^4*EA^2*alpha*cos(alpha)-4.*R^4*EA^2*sin(alpha)*cos(alpha)+4.*GIx*R^2*EA*sin(alpha)+2.*GIx*alpha*R^2*EA*cos(alpha)^2+4.*GIx*alpha*R^2*EA*cos(alpha)-4.*cos(alpha)*sin(alpha)*GIx*R^2*EA-5.*R^4*EA^2*alpha+GIx^2*alpha^3-1.*GIx^2*alpha+2.*GIx*alpha^3*R^2*EA-6.*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4.*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3);
k21(6,3)=0;
k21(6,4)=0;
k21(6,5)=0;
k21(6,6)=-1.*GIx*(GIx^2*alpha^2+R^4*EA^2*alpha^2+2.*GIx*alpha^2*R^2*EA-1.*GIx^2+2.*GIx*R^2*EA-6.*R^2*EA*alpha*sin(alpha)*GIx-6.*R^4*EA^2*alpha*sin(alpha)-8.*R^4*EA^2*cos(alpha)+R^4*EA^2*cos(alpha)^2+cos(alpha)^2*GIx^2+2.*R^2*EA*cos(alpha)*GIx*alpha^2+2.*R^4*EA^2*cos(alpha)*alpha^2-2.*R^2*EA*cos(alpha)^2*GIx+7.*R^4*EA^2)/((4.*R^4*EA^2*alpha*cos(alpha)-4.*R^4*EA^2*sin(alpha)*cos(alpha)+4.*GIx*R^2*EA*sin(alpha)+2.*GIx*alpha*R^2*EA*cos(alpha)^2+4.*GIx*alpha*R^2*EA*cos(alpha)-4.*cos(alpha)*sin(alpha)*GIx*R^2*EA-5.*R^4*EA^2*alpha+GIx^2*alpha^3-1.*GIx^2*alpha+2.*GIx*alpha^3*R^2*EA-6.*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4.*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3)*R); 

%_______________________________________________________________________________________________________________

k12=k21';

%_______________________________________________________________________________________________________________

k11(1,1)=2*(-alpha*cos(alpha)*sin(alpha)*GIx-alpha*cos(alpha)*sin(alpha)*R^2*EA+GIx*alpha^2+R^2*EA*alpha^2-2*R^2*EA*cos(alpha)^2+4*R^2*EA*cos(alpha)-2*R^2*EA)*EA*GIx/((4*R^4*EA^2*alpha*cos(alpha)-4*R^4*EA^2*sin(alpha)*cos(alpha)+4*GIx*R^2*EA*sin(alpha)+2*GIx*alpha*R^2*EA*cos(alpha)^2+4*GIx*alpha*R^2*EA*cos(alpha)-4*cos(alpha)*sin(alpha)*GIx*R^2*EA-5*R^4*EA^2*alpha+GIx^2*alpha^3-GIx^2*alpha+2*GIx*alpha^3*R^2*EA-6*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3)*R);
k11(1,2)=-2*(alpha*cos(alpha)^2*GIx+alpha*R^2*EA*cos(alpha)^2-GIx*alpha-R^2*EA*alpha-2*R^2*EA*cos(alpha)*sin(alpha)+2*R^2*EA*sin(alpha))*EA*GIx/(R*(4*R^4*EA^2*alpha*cos(alpha)-4*R^4*EA^2*sin(alpha)*cos(alpha)+4*GIx*R^2*EA*sin(alpha)+2*GIx*alpha*R^2*EA*cos(alpha)^2+4*GIx*alpha*R^2*EA*cos(alpha)-4*cos(alpha)*sin(alpha)*GIx*R^2*EA-5*R^4*EA^2*alpha+GIx^2*alpha^3-GIx^2*alpha+2*GIx*alpha^3*R^2*EA-6*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3));
k11(1,3)=0;
k11(1,4)=0;
k11(1,5)=0;
k11(1,6)=2*(-cos(alpha)^2*GIx-3*R^2*EA*cos(alpha)^2+4*R^2*EA*cos(alpha)+GIx-R^2*EA-alpha*cos(alpha)*sin(alpha)*GIx-alpha*cos(alpha)*sin(alpha)*R^2*EA+GIx*alpha^2-GIx*alpha*sin(alpha)+R^2*EA*alpha^2-R^2*EA*alpha*sin(alpha))*EA*GIx/(4*R^4*EA^2*alpha*cos(alpha)-4*R^4*EA^2*sin(alpha)*cos(alpha)+4*GIx*R^2*EA*sin(alpha)+2*GIx*alpha*R^2*EA*cos(alpha)^2+4*GIx*alpha*R^2*EA*cos(alpha)-4*cos(alpha)*sin(alpha)*GIx*R^2*EA-5*R^4*EA^2*alpha+GIx^2*alpha^3-GIx^2*alpha+2*GIx*alpha^3*R^2*EA-6*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3); 

k11(2,1)=-2*(alpha*cos(alpha)^2*GIx+alpha*R^2*EA*cos(alpha)^2-GIx*alpha-R^2*EA*alpha-2*R^2*EA*cos(alpha)*sin(alpha)+2*R^2*EA*sin(alpha))*EA*GIx/(R*(4*R^4*EA^2*alpha*cos(alpha)-4*R^4*EA^2*sin(alpha)*cos(alpha)+4*GIx*R^2*EA*sin(alpha)+2*GIx*alpha*R^2*EA*cos(alpha)^2+4*GIx*alpha*R^2*EA*cos(alpha)-4*cos(alpha)*sin(alpha)*GIx*R^2*EA-5*R^4*EA^2*alpha+GIx^2*alpha^3-GIx^2*alpha+2*GIx*alpha^3*R^2*EA-6*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3));
k11(2,2)=2*(alpha*cos(alpha)*sin(alpha)*GIx+GIx*alpha^2+R^2*EA*alpha^2+alpha*cos(alpha)*sin(alpha)*R^2*EA-2*R^2*EA+2*R^2*EA*cos(alpha)^2)*EA*GIx/(R*(4*R^4*EA^2*alpha*cos(alpha)-4*R^4*EA^2*sin(alpha)*cos(alpha)+4*GIx*R^2*EA*sin(alpha)+2*GIx*alpha*R^2*EA*cos(alpha)^2+4*GIx*alpha*R^2*EA*cos(alpha)-4*cos(alpha)*sin(alpha)*GIx*R^2*EA-5*R^4*EA^2*alpha+GIx^2*alpha^3-GIx^2*alpha+2*GIx*alpha^3*R^2*EA-6*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3));
k11(2,3)=0;
k11(2,4)=0;
k11(2,5)=0;
k11(2,6)=-2*(-GIx*cos(alpha)*sin(alpha)+alpha*cos(alpha)*GIx-2*GIx*alpha+alpha*cos(alpha)*R^2*EA-2*R^2*EA*alpha-3*R^2*EA*cos(alpha)*sin(alpha)+3*R^2*EA*sin(alpha)+alpha*cos(alpha)^2*GIx+alpha*R^2*EA*cos(alpha)^2+sin(alpha)*GIx)*EA*GIx/(4*R^4*EA^2*alpha*cos(alpha)-4*R^4*EA^2*sin(alpha)*cos(alpha)+4*GIx*R^2*EA*sin(alpha)+2*GIx*alpha*R^2*EA*cos(alpha)^2+4*GIx*alpha*R^2*EA*cos(alpha)-4*cos(alpha)*sin(alpha)*GIx*R^2*EA-5*R^4*EA^2*alpha+GIx^2*alpha^3-GIx^2*alpha+2*GIx*alpha^3*R^2*EA-6*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3);


k11(3,1)=0;
k11(3,2)=0;
k11(3,3)=alpha*EI/(R^3*(alpha^2-2+2*cos(alpha)));
k11(3,4)=(alpha-sin(alpha))*EI/(R^2*(alpha^2-2+2*cos(alpha)));
k11(3,5)=-(cos(alpha)-1)*EI/(R^2*(alpha^2-2+2*cos(alpha)));
k11(3,6)=0;

k11(4,1)=0;
k11(4,2)=0;
k11(4,3)=(alpha-sin(alpha))*EI/(R^2*(alpha^2-2+2*cos(alpha)));
k11(4,4)=(2*alpha^2-2*alpha*sin(alpha)-cos(alpha)^2+2*cos(alpha)-1)*EI/(R*alpha*(alpha^2-2+2*cos(alpha)));
k11(4,5)=-(cos(alpha)-1)*(alpha-sin(alpha))*EI/(R*alpha*(alpha^2-2+2*cos(alpha)));
k11(4,6)=0;


k11(5,1)=0;
k11(5,2)=0;
k11(5,3)=-(cos(alpha)-1)*EI/(R^2*(alpha^2-2+2*cos(alpha)));
k11(5,4)=-(cos(alpha)-1)*(alpha-sin(alpha))*EI/(R*alpha*(alpha^2-2+2*cos(alpha)));
k11(5,5)=(alpha-sin(alpha))*(alpha+sin(alpha))*EI/(R*alpha*(alpha^2-2+2*cos(alpha)));
k11(5,6)=0;

k11(6,1)=2*(-cos(alpha)^2*GIx-3*R^2*EA*cos(alpha)^2+4*R^2*EA*cos(alpha)+GIx-R^2*EA-alpha*cos(alpha)*sin(alpha)*GIx-alpha*cos(alpha)*sin(alpha)*R^2*EA+GIx*alpha^2-GIx*alpha*sin(alpha)+R^2*EA*alpha^2-R^2*EA*alpha*sin(alpha))*EA*GIx/(4*R^4*EA^2*alpha*cos(alpha)-4*R^4*EA^2*sin(alpha)*cos(alpha)+4*GIx*R^2*EA*sin(alpha)+2*GIx*alpha*R^2*EA*cos(alpha)^2+4*GIx*alpha*R^2*EA*cos(alpha)-4*cos(alpha)*sin(alpha)*GIx*R^2*EA-5*R^4*EA^2*alpha+GIx^2*alpha^3-GIx^2*alpha+2*GIx*alpha^3*R^2*EA-6*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3);
k11(6,2)=-2*(-GIx*cos(alpha)*sin(alpha)+alpha*cos(alpha)*GIx-2*GIx*alpha+alpha*cos(alpha)*R^2*EA-2*R^2*EA*alpha-3*R^2*EA*cos(alpha)*sin(alpha)+3*R^2*EA*sin(alpha)+alpha*cos(alpha)^2*GIx+alpha*R^2*EA*cos(alpha)^2+sin(alpha)*GIx)*EA*GIx/(4*R^4*EA^2*alpha*cos(alpha)-4*R^4*EA^2*sin(alpha)*cos(alpha)+4*GIx*R^2*EA*sin(alpha)+2*GIx*alpha*R^2*EA*cos(alpha)^2+4*GIx*alpha*R^2*EA*cos(alpha)-4*cos(alpha)*sin(alpha)*GIx*R^2*EA-5*R^4*EA^2*alpha+GIx^2*alpha^3-GIx^2*alpha+2*GIx*alpha^3*R^2*EA-6*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3);
k11(6,3)=0;
k11(6,4)=0;
k11(6,5)=0;
k11(6,6)=GIx*(8*R^4*EA^2*cos(alpha)+cos(alpha)^2*GIx^2+GIx^2*alpha^2+3*R^4*EA^2*alpha^2+4*GIx*alpha^2*R^2*EA-GIx^2+2*GIx*R^2*EA-4*R^4*EA^2*alpha*sin(alpha)-2*cos(alpha)*sin(alpha)*R^2*EA*alpha*GIx-2*cos(alpha)*sin(alpha)*R^4*EA^2*alpha-4*R^2*EA*alpha*sin(alpha)*GIx-R^4*EA^2-2*GIx*R^2*EA*cos(alpha)^2-7*R^4*EA^2*cos(alpha)^2)/(R*(4*R^4*EA^2*alpha*cos(alpha)-4*R^4*EA^2*sin(alpha)*cos(alpha)+4*GIx*R^2*EA*sin(alpha)+2*GIx*alpha*R^2*EA*cos(alpha)^2+4*GIx*alpha*R^2*EA*cos(alpha)-4*cos(alpha)*sin(alpha)*GIx*R^2*EA-5*R^4*EA^2*alpha+GIx^2*alpha^3-GIx^2*alpha+2*GIx*alpha^3*R^2*EA-6*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3));

%_______________________________________________________________________________________________________________

k22(1,1)=2*(-alpha*cos(alpha)*sin(alpha)*GIx-alpha*cos(alpha)*sin(alpha)*R^2*EA+GIx*alpha^2+R^2*EA*alpha^2-2*R^2*EA*cos(alpha)^2+4*R^2*EA*cos(alpha)-2*R^2*EA)*EA*GIx/(R*(4*R^4*EA^2*alpha*cos(alpha)-4*R^4*EA^2*sin(alpha)*cos(alpha)+4*GIx*R^2*EA*sin(alpha)+2*GIx*alpha*R^2*EA*cos(alpha)^2+4*GIx*alpha*R^2*EA*cos(alpha)-4*cos(alpha)*sin(alpha)*GIx*R^2*EA-5*R^4*EA^2*alpha+GIx^2*alpha^3-GIx^2*alpha+2*GIx*alpha^3*R^2*EA-6*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3));
k22(1,2)=2*(alpha*cos(alpha)^2*GIx+alpha*R^2*EA*cos(alpha)^2-GIx*alpha-R^2*EA*alpha-2*R^2*EA*cos(alpha)*sin(alpha)+2*R^2*EA*sin(alpha))*EA*GIx/(R*(4*R^4*EA^2*alpha*cos(alpha)-4*R^4*EA^2*sin(alpha)*cos(alpha)+4*GIx*R^2*EA*sin(alpha)+2*GIx*alpha*R^2*EA*cos(alpha)^2+4*GIx*alpha*R^2*EA*cos(alpha)-4*cos(alpha)*sin(alpha)*GIx*R^2*EA-5*R^4*EA^2*alpha+GIx^2*alpha^3-GIx^2*alpha+2*GIx*alpha^3*R^2*EA-6*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3));
k22(1,3)=0;
k22(1,4)=0;
k22(1,5)=0;
k22(1,6)=2*(-cos(alpha)^2*GIx-3*R^2*EA*cos(alpha)^2+4*R^2*EA*cos(alpha)+GIx-R^2*EA-alpha*cos(alpha)*sin(alpha)*GIx-alpha*cos(alpha)*sin(alpha)*R^2*EA+GIx*alpha^2-GIx*alpha*sin(alpha)+R^2*EA*alpha^2-R^2*EA*alpha*sin(alpha))*EA*GIx/(4*R^4*EA^2*alpha*cos(alpha)-4*R^4*EA^2*sin(alpha)*cos(alpha)+4*GIx*R^2*EA*sin(alpha)+2*GIx*alpha*R^2*EA*cos(alpha)^2+4*GIx*alpha*R^2*EA*cos(alpha)-4*cos(alpha)*sin(alpha)*GIx*R^2*EA-5*R^4*EA^2*alpha+GIx^2*alpha^3-GIx^2*alpha+2*GIx*alpha^3*R^2*EA-6*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3); 

k22(2,1)=2*(alpha*cos(alpha)^2*GIx+alpha*R^2*EA*cos(alpha)^2-GIx*alpha-R^2*EA*alpha-2*R^2*EA*cos(alpha)*sin(alpha)+2*R^2*EA*sin(alpha))*EA*GIx/(R*(4*R^4*EA^2*alpha*cos(alpha)-4*R^4*EA^2*sin(alpha)*cos(alpha)+4*GIx*R^2*EA*sin(alpha)+2*GIx*alpha*R^2*EA*cos(alpha)^2+4*GIx*alpha*R^2*EA*cos(alpha)-4*cos(alpha)*sin(alpha)*GIx*R^2*EA-5*R^4*EA^2*alpha+GIx^2*alpha^3-GIx^2*alpha+2*GIx*alpha^3*R^2*EA-6*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3));
k22(2,2)=2*(alpha*cos(alpha)*sin(alpha)*GIx+GIx*alpha^2+R^2*EA*alpha^2+alpha*cos(alpha)*sin(alpha)*R^2*EA-2*R^2*EA+2*R^2*EA*cos(alpha)^2)*EA*GIx/(R*(4*R^4*EA^2*alpha*cos(alpha)-4*R^4*EA^2*sin(alpha)*cos(alpha)+4*GIx*R^2*EA*sin(alpha)+2*GIx*alpha*R^2*EA*cos(alpha)^2+4*GIx*alpha*R^2*EA*cos(alpha)-4*cos(alpha)*sin(alpha)*GIx*R^2*EA-5*R^4*EA^2*alpha+GIx^2*alpha^3-GIx^2*alpha+2*GIx*alpha^3*R^2*EA-6*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3));
k22(2,3)=0;
k22(2,4)=0;
k22(2,5)=0;
k22(2,6)=2*(-GIx*cos(alpha)*sin(alpha)+alpha*cos(alpha)*GIx-2*GIx*alpha+alpha*cos(alpha)*R^2*EA-2*R^2*EA*alpha-3*R^2*EA*cos(alpha)*sin(alpha)+3*R^2*EA*sin(alpha)+alpha*cos(alpha)^2*GIx+alpha*R^2*EA*cos(alpha)^2+sin(alpha)*GIx)*EA*GIx/(4*R^4*EA^2*alpha*cos(alpha)-4*R^4*EA^2*sin(alpha)*cos(alpha)+4*GIx*R^2*EA*sin(alpha)+2*GIx*alpha*R^2*EA*cos(alpha)^2+4*GIx*alpha*R^2*EA*cos(alpha)-4*cos(alpha)*sin(alpha)*GIx*R^2*EA-5*R^4*EA^2*alpha+GIx^2*alpha^3-GIx^2*alpha+2*GIx*alpha^3*R^2*EA-6*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3);

k22(3,1)=0;
k22(3,2)=0;
k22(3,3)=alpha*EI/(R^3*(alpha^2-2+2*cos(alpha)));
k22(3,4)=(alpha-sin(alpha))*EI/((alpha^2-2+2*cos(alpha))*R^2);
k22(3,5)=EI*(cos(alpha)-1)/((alpha^2-2+2*cos(alpha))*R^2);
k22(3,6)=0;

k22(4,1)=0;
k22(4,2)=0;
k22(4,3)=(alpha-sin(alpha))*EI/((alpha^2-2+2*cos(alpha))*R^2);
k22(4,4)=EI*(2*alpha^2-2*alpha*sin(alpha)-cos(alpha)^2+2*cos(alpha)-1)/((alpha^2-2+2*cos(alpha))*alpha*R);
k22(4,5)=EI*(alpha-sin(alpha))*(cos(alpha)-1)/((alpha^2-2+2*cos(alpha))*alpha*R);
k22(4,6)=0;

k22(5,1)=0;
k22(5,2)=0;
k22(5,3)=EI*(cos(alpha)-1)/((alpha^2-2+2*cos(alpha))*R^2);
k22(5,4)=EI*(alpha-sin(alpha))*(cos(alpha)-1)/((alpha^2-2+2*cos(alpha))*alpha*R);
k22(5,5)=EI*(alpha+sin(alpha))*(alpha-sin(alpha))/((alpha^2-2+2*cos(alpha))*alpha*R);
k22(5,6)=0;

k22(6,1)=2*(-cos(alpha)^2*GIx-3*R^2*EA*cos(alpha)^2+4*R^2*EA*cos(alpha)+GIx-R^2*EA-alpha*cos(alpha)*sin(alpha)*GIx-alpha*cos(alpha)*sin(alpha)*R^2*EA+GIx*alpha^2-GIx*alpha*sin(alpha)+R^2*EA*alpha^2-R^2*EA*alpha*sin(alpha))*EA*GIx/(4*R^4*EA^2*alpha*cos(alpha)-4*R^4*EA^2*sin(alpha)*cos(alpha)+4*GIx*R^2*EA*sin(alpha)+2*GIx*alpha*R^2*EA*cos(alpha)^2+4*GIx*alpha*R^2*EA*cos(alpha)-4*cos(alpha)*sin(alpha)*GIx*R^2*EA-5*R^4*EA^2*alpha+GIx^2*alpha^3-GIx^2*alpha+2*GIx*alpha^3*R^2*EA-6*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3);
k22(6,2)=2*(-GIx*cos(alpha)*sin(alpha)+alpha*cos(alpha)*GIx-2*GIx*alpha+alpha*cos(alpha)*R^2*EA-2*R^2*EA*alpha-3*R^2*EA*cos(alpha)*sin(alpha)+3*R^2*EA*sin(alpha)+alpha*cos(alpha)^2*GIx+alpha*R^2*EA*cos(alpha)^2+sin(alpha)*GIx)*EA*GIx/(4*R^4*EA^2*alpha*cos(alpha)-4*R^4*EA^2*sin(alpha)*cos(alpha)+4*GIx*R^2*EA*sin(alpha)+2*GIx*alpha*R^2*EA*cos(alpha)^2+4*GIx*alpha*R^2*EA*cos(alpha)-4*cos(alpha)*sin(alpha)*GIx*R^2*EA-5*R^4*EA^2*alpha+GIx^2*alpha^3-GIx^2*alpha+2*GIx*alpha^3*R^2*EA-6*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3);
k22(6,3)=0;
k22(6,4)=0;
k22(6,5)=0;
k22(6,6)=(GIx^2*alpha^2+3*R^4*EA^2*alpha^2+4*GIx*alpha^2*R^2*EA+cos(alpha)^2*GIx^2+8*R^4*EA^2*cos(alpha)-GIx^2-4*R^2*EA*alpha*sin(alpha)*GIx-2*cos(alpha)*sin(alpha)*R^2*EA*alpha*GIx-2*cos(alpha)*sin(alpha)*R^4*EA^2*alpha-4*R^4*EA^2*alpha*sin(alpha)+2*GIx*R^2*EA-2*GIx*R^2*EA*cos(alpha)^2-R^4*EA^2-7*R^4*EA^2*cos(alpha)^2)*GIx/(R*(4*R^4*EA^2*alpha*cos(alpha)-4*R^4*EA^2*sin(alpha)*cos(alpha)+4*GIx*R^2*EA*sin(alpha)+2*GIx*alpha*R^2*EA*cos(alpha)^2+4*GIx*alpha*R^2*EA*cos(alpha)-4*cos(alpha)*sin(alpha)*GIx*R^2*EA-5*R^4*EA^2*alpha+GIx^2*alpha^3-GIx^2*alpha+2*GIx*alpha^3*R^2*EA-6*GIx*alpha*R^2*EA+cos(alpha)^2*R^4*EA^2*alpha+cos(alpha)^2*GIx^2*alpha+4*R^4*EA^2*sin(alpha)+R^4*EA^2*alpha^3));

%_______________________________________________________________________________________________________________

k=[k11 k12;k21 k22];

%_______________________________________________________________________________________________________________

%R=[rot2local(0,0,pi/2)',zeros(3,3);zeros(3,3),rot2local(0,0,pi/2)];
%R2=[rot2local(0,0,pi)',zeros(3,3);zeros(3,3),rot2local(0,0,pi)];
%R*k11*R'
%R'*k11*R
%R2*k11*R2'
%R2'*k11*R2

