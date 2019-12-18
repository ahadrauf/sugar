function [Cppcorj,Cjvcj,Cppj]=sc_cppcor1

tic;
j=0;
i=sqrt(-1);
for g = 1e-6 : (20-1)*1e-6/40 :10e-6
    j=j+1;
    L=10e-6; 
    [Q,Cpp,C,Cppcor,Cjvc] = sc_thin_plate_ground_1( -L + g *i, -eps + g *i);
    Qj(j)=real(Q);
    Cppj(j) = real(Cpp);
    Cj(j) = real(C);
    Cppcorj(j) = real(Cppcor);
    Cjvcj(j) = real(Cjvc);    
    gj(j) = g;
end

figure(1);plot(gj,Cppj,'r',gj,Cppcorj,'b',gj,Cjvcj,'go');grid on;

toc

title('Capacitance vs gap spacing for a 2D section of parallel plates. LWH=20um*20um*2um');
xlabel('gap spacing [m]');
ylabel('capacitance [F]: red --> inf plate approx, green --> exact, blue --> fringe correction');

%Cppcorj
%Cjvcj