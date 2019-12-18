% 3D small-deflection linear beam model
%
% Input parameters:
%   l, w  - beam length and width
%   h - beam height (often specified in process parameters)
%   density, fluid, viscosity, Youngmodulus -
%       material parameters specified as part of the process info
%
% Nodes/variables:
%   The two end nodes are conceptually in the middle of the two
%   w-by-h end faces, which are l units apart.  Each node has
%   the usual spatial displacement variables (x, y, z, rx, ry, rz)

function [output] = MF_beam3dnl6(flag, R, param, q, t, nodes, varargin);

switch(flag)

case 'vars'

  output.dynamic = {1 {'x' 'y' 'z' 'rx' 'ry' 'rz'};
                    2 {'x' 'y' 'z' 'rx' 'ry' 'rz'}};

case 'check'

  if (~isfield(param, 'density')       | ...
      ~isfield(param, 'fluid')         | ...
      ~isfield(param, 'viscosity')     | ...
      ~isfield(param, 'Youngsmodulus'))
    output = 'Missing params typically specified in process file';
  elseif ~isfield(param, 'l')
    output = 'Missing length';
  elseif ~isfield(param, 'w')
    output = 'Missing width';
  elseif ~isfield(param, 'h') 
    output = 'Missing height (layer thickness)';
  else
    output = [];
  end


case 'M'
   
   %geometrical params
   H = param.h;               %Beam layer thickness, along z-axis.
   W = param.w;               %Beam width, along y-axis.
   L = param.l;               %Beam length, along x-axis.
   D = param.density;         %Poly-Si density.
   A = H * W;                 %Beam crossectional area, yz-plane.
   Iy = W * H^3 / 12;         %Rotatory moment of inertia about y.
   Iz =	H * W^3 / 12;         %Rotatory moment of inertia about z. Used in 2D.
   Jx =	Iy + Iz;              %Polar moment of inertia, torsional inertia.
                              %J=0.14 if h<=w, cJ=1/3 if h<<w.
   A = H * W;					%Cross sectional area.   
   a =  1 / 3;					%Various quantities...
   b = 13 / 35;
   c = 13 / 35;
   d = Jx / (3 * A);
   e = L^2 / 105;
   f = L^2 / 105;
   g = 1 / 3;
   h = 13 / 35;
   i = 13 / 35;
   j = Jx / (3 * A);
   k = L^2 / 105;
   l = L^2 / 105;
   m = -11 * L / 210;
   n = 11 * L / 210;
   oo= 	9 / 70;
   p = 13 * L / 420;
   q = 9 / 70;
   r = -13 * L / 420;
   s = Jx / (6 * A);
   t = -13 * L / 420;
   u = -L^2 / 140;
   v = 1 / 6;
   w = -L^2 / 140;
   
%rotation matrix
   Rall(1:3,1:3) = R;
   Rall(4:6,4:6) = R;
   R = Rall;
  	 
%mass matrix 
  
   M11 = R * [ ...
     a     0     0     0     0     0;
     0     b     0     0     0     n;
     0     0     c     0     m     0;
     0     0     0     d     0     0;
     0     0     m     0     e     0;     
     0     n     0     0     0     f] * R';
   M22 = R * [ ...
     g     0     0     0     0     0;
     0     h     0     0     0    -n;
     0     0     i     0    -m     0;
     0     0     0     j     0     0;
     0     0    -m     0     k     0;
     0    -n     0     0     0     l] * R';
   M21 = R * [ ...
     v     0     0     0     0     0;
     0     oo    0     0     0     p;
     0     0     q     0     r     0;
     0     0     0     s     0     0;
     0     0    -r     0     w     0;
     0     t     0     0     0     u] * R';
  
   output = D*A*L*[M11  M21'; M21 M22];
  
case 'D'

%matrix params
   H = param.h;		%Beam layer thickness, along z-axis.
   W = param.w;		%Beam width, along y-axis.
   L = param.l;		%Beam length, along x-axis.
   A = H * W;		%Beam crossectional area, yz-plane.
   Mu = param.viscosity;%Viscocity of fluid environment.
   delta = param.fluid;	%Fluid layer thickness.
   Iy =	W * H^3 / 12;	%Rotatory moment of inertia about y.
   Iz =	H * W^3 / 12;	%Rotatory moment of inertia about z. Used in 2D.
   Jx =	Iy + Iz;	%Polar moment of inertia, torsional inertia.
                        %J=0.14 if h<=w, cJ=1/3 if h<<w.
   A = H * W;					%Cross sectional area.   
   a = 1 / 3;					%Various quantities...
   b = 13 / 35;
   c = 13 / 35;
   d = Jx / (3 * A);
   e = L^2 / 105;
   f = L^2 / 105;
   g = 1 / 3;
   h = 13 / 35;
   i = 13 / 35;
   j = Jx / (3 * A);
   k = L^2 / 105;
   l = L^2 / 105;
   m = -11 * L / 210;
   n = 11 * L / 210;
   oo = 9 / 70;
   p = 13 * L / 420;
   q = 9 / 70;
   r = -13 * L / 420;
   s = Jx / (6 * A);
   t = -13 * L / 420;
   u = -L^2 / 140;
   v = 1 / 6;
   w = -L^2 / 140;
   
%rotation matrix
   Rall(1:3,1:3) = R;
   Rall(4:6,4:6) = R;
   R = Rall;
   
%damping matrix. This viscous model is good for planar motion.  
   D11 = R * [ ...
     a     0     0     0     0     0;
     0     b     0     0     0     n;
     0     0     c     0     m     0;
     0     0     0     d     0     0;
     0     0     m     0     e     0;     
     0     n     0     0     0     f] * R';
   D22 = R * [ ...
     g     0     0     0     0     0;
     0     h     0     0     0    -n;
     0     0     i     0    -m     0;
     0     0     0     j     0     0;
     0     0    -m     0     k     0;
     0    -n     0     0     0     l] * R';
   D21 = R * [ ...
     v     0     0     0     0     0;
     0     oo    0     0     0     p;
     0     0     q     0     r     0;
     0     0     0     s     0     0;
     0     0    -r     0     w     0;
     0     t     0     0     0     u] * R';
  
  output = Mu*W*L/delta*[D11  D21'; D21 D22];

case 'K'

%matrix params
   E = param.Youngsmodulus;	%Modulus of elasticity.
   H = param.h;			%Beam layer thickness, along z-axis.
   W = param.w; 		%Beam width, along y-axis.
   L = param.l;			%Beam length, along x-axis.
   A = H * W;			%Beam crossectional area, yz-plane.

   Nu = param.Poisson;		%Poisson's ratio.
   G = E / (2 * (1 + Nu));	%Shear modulus. GJ = Torsional stiffness 

   if (W>H)
     J = W*H^3*(16/3-3.36*H/W*(1-(H/W)^4/12))/16;
   else
     J = H*W^3*(16/3-3.36*W/H*(1-(W/H)^4/12))/16;
   end
 
% The newer version used this different formula for J.  If the two
% are supposed to be the same, there is a very high relative error
% (at least ~0.9 for the mirror mode demo).  I don't know which one
% is correct.

%   J2 = 2/7*(W*H)^3/(W^2+H^2);	%Polar second moment of area.
%   abs((J2-J)/J)

   Iy =	W * H^3 / 12;		%Moment about y. 
   Iz =	H * W^3 / 12;		%Moment about z. Used in 2D. 

   a = E * A / L;
   b = 12 * E * Iz / L^3;
   c = 12 * E * Iy / L^3;
   d = G * J / L;
   e = 4 * E * Iy / L; %
   f = 2 * E * Iy / L;
   g = 4 * E * Iz / L;
   h = 2 * E * Iz / L;
   i = 6 * E * Iy / L^2; %
   j = 6 * E * Iz / L^2;
   
%rotation matrix
   Rall(1:3,1:3) = R;
   Rall(4:6,4:6) = R;
   R = Rall;
   
%stiffness matrix, linear model
   K11 = R * [ ...
     a     0     0     0     0     0; 
     0     b     0     0     0     j; 
     0     0     c     0    -i     0; 
     0     0     0     d     0     0; 
     0     0    -i     0     e     0;     
     0     j     0     0     0     g] * R';
   K22 = R * [ ...
     a     0     0     0     0     0; 
     0     b     0     0     0    -j; 
     0     0     c     0     i     0; 
     0     0     0     d     0     0; 
     0     0     i     0     e     0; 
     0    -j     0     0     0     g] * R';   
   K12 = R * [ ...         
    -a     0     0     0     0     0; 
     0    -b     0     0     0     j; 
     0     0    -c     0    -i     0; 
     0     0     0    -d     0     0; 
     0     0     i     0     f     0;     
     0     -j    0     0     0     h] * R';
   K21 = R * [ ...         
    -a     0     0     0     0     0; 
     0    -b     0     0     0    -j; 
     0     0    -c     0    -i     0; 
     0     0     0    -d     0     0; 
     0     0    -i     0     f     0;     
     0     j     0     0     0     h] * R';
      
  output = [K11 K21'; K21 K22];

case 'pos'

  % Compute relative positions of beam nodes

  output = R * [0 param.l;
                0 0;
                0 0];            

case 'F'

% Need dF to go with it.

   w = param.w; %width
   h = param.h; %layer thickness
   o = zeros(3,3);
   R = [R o; o R]; %rotation
   A = w*h; %cross-sectional area
   I = h^3*w/12; %Moment of inertial.
   E = param.Youngsmodulus; %Young's modulus.
   L = param.l;
   %psi, x, y, f
%[0.2, 0.98934539904854, 0.13285193550394, 0.40594289791337] 
%[0.4, 0.95752017735376, 0.26284328878478, 0.84947486001005]
%[0.6, 0.90488892146335, 0.38725775983641, 1.37915474106551]
%[0.8, 0.83185312268779, 0.50369877873648, 2.07334456831982]

ndX=[0
   0.00000008181995
   0.00011712612414
   0.00041952810875
   0.00093690044863
   0.00167076760239
   0.00259196105495
   0.00371619293222
   0.00504433635344
   0.00656111051523
   0.00829955353001
   0.01022797778489
   0.01236358524871
   0.01470744208206
   0.01724297807341
   0.01997040012606
   0.02290828848043
   0.02603921636532
   0.02938253929019
   0.03292024814157
   0.03665272646419
   0.04058039637349
   0.04472434091658
   0.04904426370881
   0.05358242267306
   0.05829696179158
   0.06320963312940
   0.06832123370147
   0.07363265406759
   0.07914489016012
   0.08483461166794
   0.09072637978448
   0.09682152067539
   0.10309530956350
   0.10957434731821
   0.11623289887706
   0.12309909737871
   0.13014614712080
   0.13740390299353
   0.14484445844364
   0.15246823163062
   0.16030803671116
   0.16830091393957
   0.17651328370428
   0.18491400825245
   0.19350443919428
   0.20232374673628
   0.21129988244238
   0.22047141119480
   0.22988248383025
   0.23945429555063
   0.24923043385576
   0.25926143861964
   0.26946055705052
   0.27987790969799
   0.29046675233087
   0.30133752719461
   0.31238987030641
   0.32368562419576
   0.33523530530740
   0.34705132441264
   0.35907643420980
   0.37139113337483
   0.38401982537642
   0.39688770681987
   0.41010202355543
   0.42369425045168
   0.43759247745833
   0.45196663104550
   0.46670208040318
   0.48189783526981
   0.49761363290914
   0.51402476633897
   0.53127334758263
   0.54929349865109
   0.56853763612495
   0.58900267342612
   0.61098034468258
   0.63562238593590
   0.65847580090150
   0.70019073653589
   0.76090960002528];
ndY=[0
   0.00066597013341
   0.01367448878652
   0.02664576251138
   0.03963575298104
   0.05264386747655
   0.06561196970744
   0.07856776791480
   0.09151040651957
   0.10440992812796
   0.11732363289344
   0.13019236250356
   0.14304445323812
   0.15587915073602
   0.16866599091088
   0.18140387873620
   0.19412174859801
   0.20678880767369
   0.21943432932652
   0.23202730071262
   0.24456674927869
   0.25705173054338
   0.26951227531504
   0.28188578797848
   0.29423349954671
   0.30649220361861
   0.31869226249492
   0.33083297098960
   0.34291367826590
   0.35493379310736
   0.36686028414849
   0.37872474416383
   0.39052676998256
   0.40223282113490
   0.41387537509252
   0.42542052319871
   0.43690140364609
   0.44828372206013
   0.45960135465772
   0.47081960384488
   0.48193786977333
   0.49299128977261
   0.50390845520695
   0.51476075252846
   0.52551230575855
   0.53616301562315
   0.54675070268997
   0.55720043487740
   0.56754984452959
   0.57783883967914
   0.58798985607236
   0.59804263423938
   0.60803974669320
   0.61790051987702
   0.62766742193404
   0.63729861339644
   0.64688360226892
   0.65633635228698
   0.66570420921784
   0.67499122232530
   0.68420225296031
   0.69329143033093
   0.70231384374923
   0.71127968339309
   0.72013538334493
   0.72894740145786
   0.73772866642772
   0.74642760882409
   0.75513737348410
   0.76378742564309
   0.77243865465610
   0.78109827713251
   0.78985886701464
   0.79877939352331
   0.80781117095708
   0.81716296719279
   0.82681190717246
   0.83687849699133
   0.84786055928501
   0.85791267152588
   0.87557417923060
   0.90093358233326];
ndpsi0=[0
   0.00063661977237
   0.01304908455367
   0.02546154933497
   0.03787401411627
   0.05028647889757
   0.06269894367886
   0.07511140846016
   0.08752387324146
   0.09993633802276
   0.11234880280406
   0.12476126758536
   0.13717373236666
   0.14958619714796
   0.16199866192926
   0.17441112671056
   0.18682359149186
   0.19923605627316
   0.21164852105446
   0.22406098583576
   0.23647345061706
   0.24888591539836
   0.26129838017966
   0.27371084496096
   0.28612330974226
   0.29853577452355
   0.31094823930485
   0.32336070408615
   0.33577316886745
   0.34818563364875
   0.36059809843005
   0.37301056321135
   0.38542302799265
   0.39783549277395
   0.41024795755525
   0.42266042233655
   0.43507288711785
   0.44748535189915
   0.45989781668045
   0.47231028146175
   0.48472274624305
   0.49713521102435
   0.50954767580565
   0.52196014058694
   0.53437260536824
   0.54678507014954
   0.55919753493084
   0.57160999971214
   0.58402246449344
   0.59643492927474
   0.60884739405604
   0.62125985883734
   0.63367232361864
   0.64608478839994
   0.65849725318124
   0.67090971796254
   0.68332218274384
   0.69573464752514
   0.70814711230644
   0.72055957708774
   0.73297204186904
   0.74538450665034
   0.75779697143163
   0.77020943621293
   0.78262190099423
   0.79503436577553
   0.80744683055683
   0.81985929533813
   0.83227176011943
   0.84468422490073
   0.85709668968203
   0.86950915446333
   0.88192161924463
   0.89433408402593
   0.90674654880723
   0.91915901358853
   0.93157147836983
   0.94398394315113
   0.95639640793243
   0.96880887271373
   0.98122133749502
   0.99363380227632];
ndFLLEI=[0
   0.00199999933939
   0.04100164175591
   0.08003563316687
   0.11913768473082
   0.15834407141507
   0.19767968958071
   0.23718191730009
   0.27688533596513
   0.31681506991213
   0.35702667847465
   0.39753540750083
   0.43839020258476
   0.47963109286027
   0.52128091232216
   0.56337717648971
   0.60598151483952
   0.64911581603389
   0.69284961957492
   0.73720467943554
   0.78222756948388
   0.82796711309791
   0.87451235004344
   0.92184483308435
   0.97010046650201
   1.01925492838833
   1.06941291265496
   1.12064229838163
   1.17301534787516
   1.22660916437686
   1.28143774686134
   1.33764767948169
   1.39533279287570
   1.45450907381340
   1.51535777805550
   1.57789703892961
   1.64233646949997
   1.70869764911867
   1.77722433710537
   1.84794338621945
   1.92099778366732
   1.99669558926351
   2.07490940884150
   2.15614255600347
   2.24044308585984
   2.32803113760823
   2.41937572862385
   2.51430530675238
   2.61332176273829
   2.71704116014021
   2.82526027946331
   2.93864648654118
   3.05802887453647
   3.18317402167578
   3.31499461338599
   3.45361142276553
   3.60072154378715
   3.75605063095309
   3.92103673774541
   4.09674664473111
   4.28442207999469
   4.48451266384277
   4.69947675968090
   4.93137802584468
   5.18093347622344
   5.45214268618877
   5.74841729064579
   6.07161085230054
   6.42932047194787
   6.82393822097680
   7.26381034910775
   7.75829603923675
   8.32316007877861
   8.97802449133790
   9.74016372471554
  10.65695623193357
  11.77167714429571
  13.16447467300540
  15.02820604064603
  17.12637943939654
  22.24082707941018
  34.98517120012139];

%Given a y displacement, what are the nonlinear restoring forces.
y=q(8)-q(2); %(+/-) net displacement from node 1 & 2. 
ndy=abs(y)/L; %nondimensional y displacement
signy=sign(y);

%Find out what the ideal nonlinear x & psi are supposed to be. 
%The interpolation works like this youtput=interp1(xdata,ydata,xinput).
X=L*interp1(ndY,ndX,ndy); %[meters]
PSI=(pi/2)*interp1(ndY,ndpsi0,ndy); %[radians]

%Find the ideal x y psi nonlinear resoring forces. 
FY=(E*I)/(L*L)*interp1(ndY,ndFLLEI,ndy); %[Newtons]
FY_minuslinear=FY-y*3*E*I/(L*L*L); %total nonlinear y-force minus the linear part.
FM_minuslinear=PSI*4*E*I/L-q(12)*4*E*I/L;
FX=X*E*A/L; 

output = [0 0 0, 0 0 0, -FX -FY_minuslinear*signy 0, 0 0 -FM_minuslinear*signy]';
   
case 'display'

  displaybeam(q, nodes(1).pos, R, param.l, param.w, param.h);
  
otherwise

  output = [];

end

