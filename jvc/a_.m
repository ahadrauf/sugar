x =[ 0     1     2
     0     1     2
     0     1     2
     0     1     2
     0     1     2];
y = [0     0     0
     1     1     1
     1     1     1
     0     0     0
     0     0     0];
z = [0     0     0
     0     0     0
     1     1     1
     1     1     1
     0     0     0];
z=[[0;0;1;1;0],[0;0;1;1;0],[0;0;1;1;0]]
y=[[0;1;1;0;0],[0;1;1;0;0],[0;1;1;0;0]]
figure(1);clf;axis equal;hold on;surfl(x+3,y,z);mesh(x+6,y,z);
figure(1);clf;axis equal;hold on;surfl(x+3,y,z);mesh(x+3,y+2,z);
surfaces={{x+3,y,z};{x+3,y+2,z}};
figure(1);rotate3d;clf;axis equal;hold on;surfl(x+3,y,z);mesh(x+3,y+2,z);
electro3d(surfaces,v)

[x,y,z]=sphere;
x=x';y=y';z=z';
z(:,1)=[];
z(:,20)=[];
x(:,1)=[];
x(:,20)=[];
y(:,1)=[];
y(:,20)=[];
figure(3);rotate3d;mesh(x,y,z);
figure(3);clf;rotate3d on;mesh(x,y,z);hold on;mesh(x+2.2,y,z);surfaces={{x,y,z},{x+2.2,y,z}};electro3d(surfaces,[10,-10]);

