%Find t given z for the Schwarz-Christofel transformation
%of a comb drive corner
function [result]=sc_comb1(t,p,g,z)
%By: Jason Vaughn Clark - Nov2001
u=sqrt((t-g^2/(p^2))/(t+1)); %Parameter for SC
zero=-z+(2*p*atan(p*u/g)+log(-(1+u)/(-1+u))*g)/pi; %Schwarz-Christofel transformation
result=real(zero)+imag(zero); %add the real and imaginary contributions together; looking for 0+0. Later check for R-C=0!

%figure(1);clf;hold on;for t=-1+eps:-0.1:-100,plot(t,sc_comb1(t));end
%p=2e-6;g=2e-6;z=3e-6+j*g;tv=[-1+eps,-1e-10];th=[-1-eps,1e10]; [T,Zero]=fzero('sc_comb1',tv,[],p,g,z);

