i=0;
i=i+1;load(i) =...
   39.2700;
dc(i) =...
   20.3300;
dsp(i) =...
   26.9100;
dof(i) =...
        5226;
     
     i=i+1;load(i) =...
   37.8400;
dc(i) =...
   11.8700;
dsp(i) =...
   20.9800;
dof(i) =...
        4026;
 
 
 i=i+1;load(i) =...
   19.8300;
dc(i) =...
    6.2600;
dsp(i) =...
   13.5700;
dof(i) =...
        2826;
 
 
 i=i+1;load(i) =...
    8.5700;
dc(i) =...
    2.6400;
dsp(i) =...
    7.9100;
dof(i) =...
        1626;
 
 
 i=i+1;load(i) =...
    3.7900;
dc(i) =...
    1.2600;
dsp(i) =...
    4.4500;
dof(i) =...
906

i=i+1;load(i) =...
   52.2300;
dc(i) =...
    29
dsp(i) =...
    34
dof(i) =...
        6426;
     
     i=i+1;load(i) =...
  128.5800;
dc(i) =...
   63.4400;
dsp(i) =...
   43.1200;
dof(i) =...
7626;

i=i+1;load(i) =...
  194.8200;
dc(i) =...
  210.4700;
dsp(i) =...
   50.3700;
dof(i) =...
8826;

i=i+1;load(i) =...
  287.4200;
dc(i) =...
  474.7800;
dsp(i) =...
   59.2100;
dof(i) =...
10026;

DOF=dof;
ldof=length(dof);
DOF=sort(dof);
for j = 1 : ldof
   k(j) = find(DOF(j)==dof);
end
clf;
plot(dof(k),load(k)/60,'r',dof(k),dc(k)/60,'g',dof(k),dsp(k)/60,'b')
hold on;
plot(dof(k),load(k)/60,'*r',dof(k),dc(k)/60,'*g',dof(k),dsp(k)/60,'*b')
hold off;
%plot(dof(k),load(k),'r')

