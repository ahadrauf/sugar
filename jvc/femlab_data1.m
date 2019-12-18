%femlab_parallelplate3
plate_length = 20;
plate_width = 20;
plate_thickness = 2;
gap_space_between_plates = 2;
micron_scale = 1e-6;
boundarymargin = 20;

course -> fine

capacitance_femlab(1) =  2.5112e-015
capacitance_parallel_approx(1) =  1.7708e-015
elapsed_time =   25.6670

capacitance_femlab(2) =  2.5058e-015
capacitance_parallel_approx(2) =  1.7708e-015
elapsed_time =   47.0980
   
capacitance_femlab(3) =  2.5011e-015
capacitance_parallel_approx(3) =  1.7708e-015
elapsed_time =   84.1610

capacitance_femlab(4) =  2.4969e-015
capacitance_parallel_approx(4) =  1.7708e-015
elapsed_time =  320.8910
  
  
================================================
  
cap(1) =  2.1586e-015
capacitance_pp(1) =  1.7708e-015
ans(1) =    0.0500
cap(2) =  2.1582e-015
capacitance_pp(2) =  1.7708e-015
ans(2) =    0.1000
cap(3) =  2.1538e-015
capacitance_pp(3) =  1.7708e-015
ans(3) =    0.5000
cap(4) =  2.1467e-015
capacitance_pp(4) =  1.7708e-015
ans(4) =     1
cap(5) =  2.1328e-015
capacitance_pp(5) =  1.7708e-015
ans(5) =     2

======================

fringing_cap = sum(-Q{1}([lhs max(lhs)+1]) + Q{1}(640))/15  + 2.1586e-015 + sum(-Q{1}([rhs max(rhs)+1]) + Q{1}(220))/15 + 2*sum(-Q{1}([ths max(ths)]))/15
fringing_cap = 2.4634e-015

===================

>> sum( -Q{1}([lhs max(lhs)+1])*2 + Q{1}(640) )/15 / fringing_cap
ans =    0.0853
>> sum( -Q{1}([rhs max(rhs)+1])*2 + Q{1}(220) )/15 / fringing_cap
ans =    0.7474
>> sum( -Q{1}([ths max(ths)])*4 )/15 / fringing_cap
ans =    0.1703
>> 0.1703+0.7474+0.0853
ans =    1.0030

===============================================

>> fringing_cap / femlab
2.4634e-015 / 2.4969e-015
ans = 0.9866


clf;for V=0:0.5:10 ;V0=10;G=10;b=5;a=0.5;x=0:0.1:20;y=V*(G+b*sin(a*x))/V0;figure(1);plot(0,0,0,G,x,y);grid on; hold on; end
