function F = funcAC(a)
global gamma q beta w0;
F(1)=((1+3*gamma/(8*w0^2)*a(1)^2+sqrt((q/(2*w0^2*a(1)))^2-beta^2)))-((1+3*gamma/(16*w0^2)*a(1)^2+9*gamma/(16*w0^2)*a(1)^2)+sqrt(((1+3*gamma/(16*w0^2)*a(1)^2+9*gamma/(16*w0^2)*a(1)^2))^2-((1+3*gamma/(8*w0^2)*a(1)^2)*(1+9*gamma/(8*w0^2)*a(1)^2)+beta^2)));
F(2)=((1+3*gamma/(8*w0^2)*a(2)^2-sqrt((q/(2*w0^2*a(2)))^2-beta^2)))-((1+3*gamma/(16*w0^2)*a(2)^2+9*gamma/(16*w0^2)*a(2)^2)-sqrt(((1+3*gamma/(16*w0^2)*a(2)^2+9*gamma/(16*w0^2)*a(2)^2))^2-((1+3*gamma/(8*w0^2)*a(2)^2)*(1+9*gamma/(8*w0^2)*a(2)^2)+beta^2)));

