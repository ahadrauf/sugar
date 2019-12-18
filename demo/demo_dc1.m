%Sugar tutorial of static analysis of a single beam.
%by JVClark, Oct2000, 
%Modified for new netlist format by Ningning Zhou, 02/07/2001

%intro page
clc
disp(' ');
disp('SUGAR demo - Static Analysis of a Cantilever Beam');
disp(' ');
disp('This demo consists of:');
disp('  1) Netlist of a cantilever beam subjected to an external force'); 
disp('  2) Static analysis of deflection');
disp('  3) Display of structure');
disp('  4) Displacement values');
disp(' ');
disp('Press any key to continue.'); 
pause;

%1) netlist page
clc
disp(' ');
disp('1/4) Netlist file: cantilever.net');
disp(' ');
disp('Press any key to display netlist');
pause;
disp(' ');
type cantilever.net
disp(' ');
disp('Press any key to continue.');
pause;

%2) static analysis page
clc;
disp(' ');
disp('2/4) Static analysis for deflection.');
disp(' ');
disp('These Matlab commands will load the netlist, then calculate its') 
disp('deflection:');
disp(' ');
disp('net = cho_load(''cantilever.net'')');
disp(' ');
disp('dq = cho_dc(net);');
disp(' ');
net = cho_load('cantilever.net');
dq = cho_dc(net);
disp('Press any key to continue.'); 
pause;

%3) display page
clc;
disp(' ');
disp('3/4) Display of structure.');
disp(' ');
disp('The command');
disp(' cho_display(net);');
disp('displays the undeflected structure of the netlist');
disp(' ');
disp('The command');
disp(' cho_display(net,dq);');
disp('displays the undeflected structure of the netlist');
disp(' ');
figure(1); cho_display(net);    title('Undeflected structure');
figure(2); cho_display(net,dq); title('Deflected structure');
disp('Press any key to continue');
pause;

%4) Displacement value page
clc
disp(' ');
disp('4/4) Displacement value.');
disp(' ');
disp('To find the displacement of a coordinate at a particular node');
disp('(for example, the y displacement of the cantilever tip), use');
disp(' ');
disp(' dy = cho_dq_view(dq,net,''tip'',''y'')');
disp(' ');
disp('where tip is the node name, and y is desired coordinate');
disp('For this example, the computed deflection is');
disp(' ');
dy = cho_dq_view(dq,net,'tip','y')
disp(' ');
disp('Note: Though netlist angles are in degrees, output displacement');
disp('      angles are in radians. xyz displacements are meters.');
