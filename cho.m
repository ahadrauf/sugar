function cho
% SUGAR
%   Sugar - Tool for design, modeling, and simulation of Microelectro
%   mechanical systems (MEMS). MEMS can be designed in Sugar using Sugar's
%   netlist description language. A MEMS designed in Sugar is represented
%   by its netlist.
%   
%   SUGAR COMMANDS
%
%   net = cho_load('netlistname.m') loads a netlist name
%
%   net = cho_load('netlistname.m',param) loads a netlist with parameters
%   defined in param filed. See examples below
%
%   q = cho_dc(net) computes the static displacement of the in the netlist. 
%   q contains the displacement of each degree of freedom in MEMS
%
%   cho_display(net) displays the 3D structure of MEMS
%
%   cho_display(net,q) display the displaced 3D structure of MEMS. The
%   total displacement of the MEMS is colormapped on to the surface of the
%   displaced structure
%   
%   id = lookup_coord(net,node,coord) looks up the index of a degree of 
%   freedom of a given node and coordinate in the local state vector
%   
%   q(id) can be used to obtain the displacement of a given degree of
%   freedom where id is obtained from id = lookup_coord(net,node,coord)
%
%   shownodes(net) displays the node names of the MEMS at their respective 
%   locations
%
%   TOOLS RELATED TO SUGAR
%   
%   SugarCube - Online tool for parameterization and optimization of
%   ready-made MEMS. The ready-made MEMS are available in the hierarchical
%   library inside SugarCube.
%
%   Optimization in SugarCube returns the properties (geometric, material,
%   and dynamic) of the MEMS as a function of performance. Currently
%   SugarCube can optimize the MEMS as a function of static analysis and
%   resonant frequency of a given eigen mode
%
%   CHO2GDSII is a tool that generates a GDSII layout file from Sugar
%   model. GDSII is a popular layout description format used by many 
%   designers across the world. These layouts are used to generate the 
%   masks required for microfabrication process. This tool is integrated in
%   SugarCube.
%
%   SugarCube is available online at www.nanoHUB.org/resources/sugarcube
%
%   SugarX is an online tool that connects experiment to simulation. In 
%   SugarX user directly controls the MEMS in laboratory from nanoHUB. User
%   provides the input voltages for actuating the device. SugarX collects
%   the required measurements at these voltages, uses them to calibrate the
%   device (Electro-micro metrology), and extracts the actual properties of
%   the device. These extracted properties are then substituted into the
%   simulation. This allows for reasonable comparison of experiment to
%   simulation.
%
%   iSugar is an framework that integrates Sugar, COMSOL, and SIMULINK.
%   Sugar is known for its ease of use and quickness in results; COMSSOL is
%   a tool for multi energy interaction finite element analysis (FEA);
%   SIMULINK is a tool for gaprhical building-block modeling of dynamic
%   systems. In iSugar, all these three tools are integrated thereby adding
%   up the advantages of each individual tool

help cho