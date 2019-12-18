I=sqrt(-1);
Z=-10e-6 + I*gap;
width=1e-6;
gap=1e-6;

%given z, find t
P=width;
G=gap;
K='(1+(P)/(G))';
A=strrep('(-1+2*(K)^2+ 2*(K)*sqrt((K)^2-1))','K',K);
R='(sqrt((t+1)/(t+(A))))';
realzoft_0=strrep('real((G)/pi*(((A)+1)/sqrt(A)*atanh(R)+((A)-1)/sqrt(A)*(R)/(1-(R)^2)+log(((R)*sqrt(A)-1)/((R)*sqrt(A)+1))) - (Z))','R',R);  
realzoft_0=strrep(realzoft_0,'A',A);
realzoft_0=strrep(realzoft_0,'Z',num2str(Z));
realzoft_0=strrep(realzoft_0,'P',num2str(P));
realzoft_0=strrep(realzoft_0,'G',num2str(G));

t=-1.240862105804430e-014;
eval(realzoft_0)
f=inline(realzoft_0);
T=fzero(f,[-0.8 -1e-100],optimset('disp','iter'))



