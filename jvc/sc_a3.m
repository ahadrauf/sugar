
f = rectmap(drawpoly);
r = rectangle(f);
M = pi/max(imag(r));
a = exp(min(real(r))*M);
b = exp(max(real(r))*M);
rad = log(linspace(a,b,8))/M;
plot(f,rad(2:end-1),6);
H = copyobj(get(gca,'child'),gca);
for n=1:length(H),
    set(H(n),'xdata',-get(H(n),'xdata'))
end
set(H,'linesty',' - -'), axis auto
