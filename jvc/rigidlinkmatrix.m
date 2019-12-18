function T=rigidlinkmatrix(params,R)
%By Jason Vaughn Clark, Nov/2001
if isfield(params,'L1'),L1=params.L1;else,L1=0;end
if isfield(params,'ox1'),ox1=params.ox1;else,ox1=0;end
if isfield(params,'oy1'),oy1=params.oy1;else,oy1=0;end
if isfield(params,'oz1'),oz1=params.oz1;else,oz1=0;end
if isfield(params,'L2'),L2=params.L2;else,L2=0;end
if isfield(params,'ox2'),ox2=params.ox2;else,ox2=0;end
if isfield(params,'oy2'),oy2=params.oy2;else,oy2=0;end
if isfield(params,'oz2'),oz2=params.oz2;else,oz2=0;end
r1=rot2global(ox1,oy1,oz1);
r2=rot2global(ox2,oy2,oz2);
t1=R*r1*-[0,0,0;0,0,L1;0,-L1,0]*r1'*R';
t2=R*r2* [0,0,0;0,0,L2;0,-L2,0]*r2'*R';
i3=eye(3);
o3=zeros(3);
o6=zeros(6);
tau1=[i3,o3;t1,i3];
tau2=[i3,o3;t2,i3];
T=[tau1,o6;o6,tau2];

