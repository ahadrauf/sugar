function harm_undcoef
b=8/3;
a=1;
k=0.2;
B=1;
i=0;
I=sqrt(-1);
for w=0:0.01:4
%   (A1/B)^2=1/((a-w^2+0.75*b*A1^2)^2+(k*w)^2); 
    i=i+1;
    W(i)=w;
    A1(i)=-1/3*sqrt(2)*sqrt(b*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3)*(-12*k^2*w^2+(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(2/3)+4*a^2-8*a*w^2+4*w^4+4*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3)*w^2-4*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3)*a))/(b*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3));
    A2(i)=-1/3*sqrt(b*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3)*(-(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(2/3)-4*a^2+8*a*w^2+12*k^2*w^2-4*w^4+8*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3)*w^2-8*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3)*a-I*sqrt(3)*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(2/3)+4*I*sqrt(3)*a^2-8*I*sqrt(3)*a*w^2-12*I*sqrt(3)*k^2*w^2+4*I*sqrt(3)*w^4))/(b*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3));
    A3(i)=-1/3*sqrt(-b*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3)*((-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(2/3)+4*a^2-8*a*w^2-12*k^2*w^2+4*w^4-8*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3)*w^2+8*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3)*a-I*sqrt(3)*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(2/3)+4*I*sqrt(3)*a^2-8*I*sqrt(3)*a*w^2-12*I*sqrt(3)*k^2*w^2+4*I*sqrt(3)*w^4))/(b*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3));
    A4(i)=1/3*sqrt(2)*sqrt(b*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3)*(-12*k^2*w^2+(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(2/3)+4*a^2-8*a*w^2+4*w^4+4*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3)*w^2-4*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3)*a))/(b*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3));    
    A5(i)=1/3*sqrt(-b*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3)*((-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(2/3)+4*a^2-8*a*w^2-12*k^2*w^2+4*w^4-8*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3)*w^2+8*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3)*a-I*sqrt(3)*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(2/3)+4*I*sqrt(3)*a^2-8*I*sqrt(3)*a*w^2-12*I*sqrt(3)*k^2*w^2+4*I*sqrt(3)*w^4))/(b*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3));
    A6(i)=1/3*sqrt(b*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3)*(-(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(2/3)-4*a^2+8*a*w^2+12*k^2*w^2-4*w^4+8*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3)*w^2-8*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3)*a-I*sqrt(3)*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(2/3)+4*I*sqrt(3)*a^2-8*I*sqrt(3)*a*w^2-12*I*sqrt(3)*k^2*w^2+4*I*sqrt(3)*w^4))/(b*(-24*w^2*a^2+8*a^3+24*w^4*a-72*k^2*w^4+72*k^2*w^2*a-8*w^6+81*B^2*b+3*sqrt(3)*sqrt(64*a^4*k^2*w^2+64*k^6*w^6+128*k^4*w^8+64*k^2*w^10-256*a^3*w^4*k^2+128*a^2*k^4*w^4+384*a^2*k^2*w^6-256*a*w^6*k^4-256*a*w^8*k^2+48*a^3*B^2*b-48*w^6*B^2*b+144*w^4*a*B^2*b-432*k^2*w^4*B^2*b+243*B^4*b^2-144*w^2*a^2*B^2*b+432*k^2*w^2*a*B^2*b))^(1/3));

end
subplot(1,1,1);
plot(W,A1);
subplot(1,2,1);
plot(W,A2);
 figure(3);
plot(W,A3);
 figure(4);
plot(W,A4);
 figure(5);
plot(W,A5);
 figure(6);
plot(W,A6);
