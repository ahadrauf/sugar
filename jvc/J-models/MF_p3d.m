function [output] = MF_p3d(flag,R,param,x,t,nodes,varargin)
switch(flag)
case 'vars'
  output.dynamic = {1 {'x' 'y' 'z' 'rx' 'ry' 'rz'}; 2 {'x' 'y' 'z' 'rx' 'ry' 'rz'}};
case 'F'
   if isfield(param, 'P') %Load 
      P=param.P;
      L=param.l;
      psi1=strrep('.*(1-3.*(x./L).^2+2.*(x./L).^3)','L',num2str(L));
      psi2=strrep('.*(x.*(1-x./L).^2)','L',num2str(L));
      psi3=strrep('.*(3.*(x./L).^2-2.*(x./L).^3)','L',num2str(L));
      psi4=strrep('.*(x.^2./L.*(x./L-1))','L',num2str(L));
      P1=quad8(inline(strcat(P,psi1)),0,L);
      P2=quad8(inline(strcat(P,psi2)),0,L);
      P3=quad8(inline(strcat(P,psi3)),0,L);
      P4=quad8(inline(strcat(P,psi4)),0,L);
      r=eye(3); o=zeros(3,3);
      R=[r o; o r];
      F1=R*[0;P1;0; 0;0;P2];
      F2=R*[0;P3;0; 0;0;P4];
      output=[F1;F2];
   else
      output=zeros(6,1);
   end
otherwise
   output = [];
end

