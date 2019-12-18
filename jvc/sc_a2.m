clf;
p=polygon([i,-i,Inf],[3/2,1/2,-1]);
p=polygon([2i,  -i,   1-i, 1-2i, 2-2i, 2-i,    Inf, ],... % 3i-1, 4i-1, Inf], ...
          [3/2, 1/2,  3/2, 1/2,  1/2,  3/2,    -1]);%,    3/2,  3/2,  -1]);
f=hplmap(p);
axis([-3 3 -1.5 4.5]), hold on
plot(f, 0.7*(-10:6),0.7*(1:12))

%p=polygon([-5-i,-5-3i,5-3i,5+i,4+3i,-5+3i]);
%p=polygon([-Inf+i, i, -i, Inf-i, Inf+3i, -Inf+3i]);
%f=rectmap(p,[1 2 3 4]);
%plot(f);
