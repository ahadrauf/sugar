function harm_nonlin1

lambda=5;
A=1;
w0=1;
wf=5;
i=0;
for w=0:wf/100:wf
    i=i+1;
    W(i)=w;%0*(1-3*lambda*A^2/(4*w0^2))^0.5;
    A(i)=((-w^2+w0^2)*4*w0^2/3/lambda)^0.5;
end
figure(1);plot(W,A);grid on;
