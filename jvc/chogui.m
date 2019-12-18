function chogui
clf; grid on; hold on; rotate3d;
maxaxis=500e-6;
axis([0 maxaxis 0 maxaxis]);
button=1; i=0;
while(button~=3)
   [x,y,button]=ginput(2); i=i+2;
   parameter.w=20e-6; parameter.h=5e-6;
   parameter.l=sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);
   Rp=[x(1);y(1);0];
   parameter.oz=asin((y(2)-y(1))/parameter.l);
   parameter.ox=0;parameter.oy=0;
   q_global=zeros(12,1);
   R=rot2local(0,0,parameter.oz)';
   displaybeam(q_global, Rp, R, parameter.l, parameter.w, parameter.h);
   fprintf('beam3d p1 [%d %d][l=%f w=%f oz=%f]\n',i-1,i,parameter.l,parameter.w,parameter.oz);
end
