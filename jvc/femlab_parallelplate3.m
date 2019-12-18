%femlab_parallelplate3
%Femlab code for the capacitance of 2 parallel plates
%Programmed by JVCLark - summer 2003

%Define plate geometry
plate_length = 20;
plate_width = 20;
plate_thickness = 2;
gap_space_between_plates = 2;
micron_scale = 1e-6;
boundarymargin = 20;

%Start the modeling by clearing the fem variable in the Matlab workspace.
clear fem

%Begin by specifying space dimensions and defining a geometry.
fem.sdim = {'x','y','z'};


g1=block3(plate_length,plate_width,plate_thickness,'corner',[0 0  0],[0 0 1],0);
g2=block3(plate_length,plate_width,plate_thickness,'corner',[0 0 -(gap_space_between_plates+plate_thickness)],[0 0 1],0);

co1 = geomcomp({g1}, 'face', 'all', 'edge', 'all');
co2 = geomcomp({g2}, 'face', 'all', 'edge', 'all');

g3=block3(plate_length+boundarymargin, plate_width+boundarymargin, plate_thickness*2+gap_space_between_plates+boundarymargin, ...
    'corner',[-boundarymargin/2 -boundarymargin/2 -(boundarymargin/2 + gap_space_between_plates + plate_thickness) ],[0 0 1],0);

co3 = geomcomp({co1  co2 g3}, 'ns', {'CO1', 'CO2', 'BLK1'},'sf','BLK1 - (CO1 + CO2)', 'face', 'all', 'edge', 'all');

fem.geom = scale(co3, micron_scale, micron_scale, micron_scale, 0,0,0);

%The application mode sructure contanins all the physical properties of the model. Create this with the following commands:
clear appl

oApp = flcemqv3d;
appl.mode = oApp;

appl.equ.rho0 = '0';
appl.equ.epsilonr = '1'; %'4.2';

appl.bnd.V = {'0','V0','0'};
appl.bnd.type = {'n0', 'V', 'V0'};

appl.bnd.ind = ones(1,18)*2;
appl.bnd.ind(1, [1:5 18]) = 1;
appl.bnd.ind(1, [6:9 14 16]) = 3;

fem.appl = appl;

fem = multiphysics(fem);

%Define the electric potential as a variable.
fem.variables.V0 = 15;

%Create a mesh for the geometry with the following commands. The  extended mesh structure contains all the necessary data defining the different elements.
fem.mesh = meshinit(fem, 'Hmaxfact', 1.5, 'Hgrad', 1.5, 'Hcurve', 0.5, 'Hcutoff', 0.02 ); %course
%fem.mesh = meshinit(fem, 'Hmaxfact', 1, 'Hgrad', 1.4, 'Hcurve', 0.4, 'Hcutoff', 0.01 ); %normal
%fem.mesh = meshinit(fem, 'Hmaxfact', 0.8, 'Hgrad', 1.35, 'Hcurve', 0.35, 'Hcutoff', 0.008 ); %fine
%fem.mesh = meshinit(fem, 'Hmaxfact', 0.55, 'Hgrad', 1.35, 'Hcurve', 0.3, 'Hcutoff', 0.005 ); %finer

fem.xmesh = meshextend(fem);

%Solve the problem and visualize the result.
fem.sol = femiter(fem, 'report', 'on', 'prefun', 'luinc');

bdlList = setdiff(1:18, [1  2   4]);
figure(1);
postplot(fem, 'bdl', bdlList, 'tridata',{'V','cont','on'}, 'tribar','on','arrowdata',{'Ex','Ey','Ez'},'axisequal','on','axisvisible','off');

%The capacitance can be computed using subdomain integration.

capacitance_femlab = postint(fem, '2*We/V0^2','dl',1)
capacitance_parallel_approx = 8.854187817e-12 * plate_length * plate_width / gap_space_between_plates * micron_scale 

