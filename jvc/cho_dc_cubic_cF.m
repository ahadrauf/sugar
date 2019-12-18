
function [qlin,qnl]=cho_dc_cubic_cF(net,ka,kb,F)
qlin=ka\F;
%qnl=1/6*((108*F+12*sqrt(3)*sqrt((4*ka.^3+27*kb*F.^2)*inv(kb)))*kb.^2).^(1/3)*inv(kb)-2*ka*inv(((108*F+12*sqrt(3)*sqrt((4*ka.^3+27*F.^2*kb)/kb))*kb^2).^(1/3));
1/6*108*F;
12*sqrt(3);
sqrt((4*ka.^3+27*diag(F.^2)*kb));
(inv(kb)*kb.^2).^(1/3)*inv(kb);
-2*ka*inv(((108*F+12*sqrt(3)*sqrt((4*ka.^3+27*F.^2*kb)/kb))*kb^2).^(1/3));



