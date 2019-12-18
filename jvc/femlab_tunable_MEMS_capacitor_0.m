%femlab_tunable_MEMS_capacitor_0
%This program was taken straight out of Femlab's Electromagnetics book, page210.
%June 5, 2003.

%Start the modeling by clearing the fem variable in the Matlab workspace.
clear fem

%Begin by specifying space dimensions and defining a geometry.
fem.sdim = {'x','y','z'};
blks{1} = block3(22,60,8,'corner',[0    240 46]);
blks{2} = block3(40,22,8,'corner',[22   259 46]);
blks{3} = block3(176,262,8 ,'corner',[62    19  46]);
blks{4} = block3(40,22,8 ,'corner',[238 259 46]);
blks{5} = block3(22,60,8 ,'corner',[278 240 46]);
blks{6} = block3(40,22,8 ,'corner',[238 19  46]);
blks{7} = block3(22,60,8 ,'corner',[278 0   46]);
blks{8} = block3(40,22,8 ,'corner',[22  19  46]);
blks{9} = block3(22,218,8 ,'corner',[0  41  0]);
blks{10} = block3(40,22,8 ,'corner',[-40    139  0]);
cyl1 = cylinder3(5.5,   38  ,   [11 250 8]);
co1 = geomcomp({blks{:} cyl1}, 'face', 'all', 'edge', 'all');

blk1 = block3(176,262,8,'corner',[62    19  8]);
blk2 = block3(181,22,8,'corner',[139    139 0]);
co2 = blk1 + blk2;

blk3 = block3(360,340,94,'corner',[-40  -20 -21]);
co3 = geomcomp({co1  co2 blk3}, 'ns', {'CO1', 'CO2', 'BLK1'},'sf','BLK1 - (CO1 + CO2)', 'face', 'all', 'edge', 'all');

fem.geom = scale(co3, 1e-6, 1e-6, 1e-6, 0,0,0);

%The application mode sructure contanins all the physical properties of the model. Create this with the following commands:
clear appl

oApp = flcemqv3d;
appl.mode = oApp;

appl.equ.rho0 = '0';
appl.equ.epsilonr = '4.2';

appl.bnd.V = {'0','V0','0'};
appl.bnd.type = {'n0', 'V', 'V0'};

appl.bnd.ind = ones(1,76)*2;
appl.bnd.ind(1, [1:4    9   76]) = 1;
appl.bnd.ind(1, [37:40  45  48:51   57]) = 3;

fem.appl = appl;

fem = multiphysics(fem);

%Define the electric potential as a variable.
fem.variables.V0 = 1e-6;

%Create a mesh for the geometry with the following commands. The  extended mesh structure contains all the necessary data defining the different elements.
fem.mesh = meshinit(fem, 'Hcurve', 0.5, 'Hcutoff', 0.02, 'Hgrad', 1.5, 'Hmaxfact', 1.5);

fem.xmesh = meshextend(fem);

%Solve the problem and visualize the result.
fem.sol = femiter(fem, 'report', 'on', 'prefun', 'luinc');

bdlList = setdiff(1:76, [1  2   4]);
postplot(fem, 'bdl', bdlList, 'tridata',{'V','cont','on'}, 'tribar','on','arrowdata',{'Ex','Ey','Ez'},'axisequal','on','axisvisible','off');

%The capacitance can be computed using subdomain integration.

capacitance = postint(fem, '2*We/V0^2','dl',1)

