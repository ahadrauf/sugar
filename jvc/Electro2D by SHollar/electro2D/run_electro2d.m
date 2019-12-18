
file1 = 'seg1.m';
file2 = 'seg2.m';
seg1 = load_and_parse(file1);
seg2 = load_and_parse(file2);

delta_approx = 2 ; %microns
permitivity = 8.854e-6;% mks units * 1e6 (avoids truncation error)
z_depth = 30 ; % microns

[capacitance, FF_X, FF_Y] = electro2d(seg1,seg2,...
	 delta_approx, permitivity, z_depth)

%figure(2);
%plot_electro2d(seg1, seg2, delta_approx);
