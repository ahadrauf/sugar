net=cho_load('condensation1.m');[q,k,f]=cho_dc2(net);figure(1);cho_display(net,q);

kbb=k(1:17,1:17);
kbc=k(1:17,18);
kcc=k(18,18);
kcb=k(18,1:17);
khatcc=kcc-kcb*inv(kbb)*kbc;
k2=kcb*inv(kbb);
fc=f(18);
fb=f(1:17,1);
flops(0);q2=inv(khatcc)*(fc-k2*fb);flops
flops(0);q18=inv(k)*f;flops

kbb=k(1:17,1:17);
kbc=k(1:17,18);
kcc=k(18,18);
kcb=k(18,1:17);
khatcc=sparse(kcc-kcb*inv(kbb)*kbc);
k2=sparse(kcb*inv(kbb));
fc=f(18);
fb=sparse(f(1:17,1));
K=sparse(k);
F=sparse(f);
flops(0);inv(khatcc)*(fc-k2*fb),flops
flops(0);q18=inv(K)*F;flops
