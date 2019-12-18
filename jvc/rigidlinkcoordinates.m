function [RL1,RL2,r1,r2]=rigidlinkcoordinates(params);
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
RL1=r1*[-L1;0;0];
RL2=r2*[L2;0;0];
