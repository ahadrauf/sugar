function y=j1(x)

y=5*x+1;

figure(1);hold on;grid on;

plot(x,y,'go');


for x=-2:0.01:30
    y=1*x+0;
    plot(x,y,'b');
end

