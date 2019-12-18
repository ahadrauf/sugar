%
%       This script is an extension of electro2d.
%       sub_force_matrices includes part of the calculations for
%       determining the force of one line segment on another.
%


Y_rr=-alpha*sinth+R_perp;
RTH=alpha*costh+R_par;
Y2_rr=Y_rr.^2;
X_alpha=X_trans+alpha;

LRTH=L+RTH;
LOG_value=1/2*log(LRTH.^2+Y2_rr+delta);

sub_x_plusL=LRTH.*LOG_value...
   +Y_rr.*atan(LRTH./(Y_rr+delta))-LRTH;
sub_y_plusL=LRTH.*less_zero_plusL...
   .*atan((X_alpha+L*costh)./(abs(denom_plusL)...
   +delta))-Y_rr.*LOG_value;

X_alpha+L*costh;

LRTH.*less_zero_plusL...
   .*atan((X_alpha+L*costh)./(abs(denom_plusL)...
   +delta));

y_discont=(-Y_trans.*div_sinth+RTH)*pi...
   .*sign(X_trans+alpha-Y_cos.*div_sinth);

LRTH=-L+RTH;
LOG_value=1/2*log(LRTH.^2+Y2_rr+delta);

sub_x_minusL=LRTH.*LOG_value...
   +Y_rr.*atan(LRTH./(Y_rr+delta))-LRTH;
sub_y_minusL=LRTH.*less_zero_minusL...
   .*atan((X_alpha-L*costh)./(abs(denom_minusL)...
   +delta))-Y_rr.*LOG_value;

sub_x=sub_x_plusL-sub_x_minusL;
sub_y=sub_y_plusL-sub_y_minusL;

clear Y_rr RTH Y2_rr X_alpha LRTH LOG_value ...
	sub_x_plusL sub_x_minusL sub_y_plusL ...
	sub_y_minusL


