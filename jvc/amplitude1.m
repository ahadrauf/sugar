
%natural_frequencies
%modeshapes
%normalized_modeshape_matrix
%decoupled_equations
%general_solution
%maximum_response
%quality_factor

%M*qdotdot + C*qdot + K*dot = F;
 


%[X,ww] = eig(K,M);


%m=X'*M*X;
%c=X'*C*X;
%k=X'*K*X;

%w=[18.850;42.390;58.718];
%c=2*0.02*w;
clear all;
figure(1);
hold on;
i=0;
for w = 1:100 
   i=i+1;
   c=2*0.02*[18.858,0,0;0,42.390,0;0,0,58.718];
   m=[1,0,0;0,1,0;0,0,1];
   k=[355.31,0,0;0,1796.87,0;0,0,3447.31];
   F0=[4.2336;-4.4926;3.3105] * 1e-3 * 600;
   A = ((k-m*w.^2).^2 + c.^2 * w.^2).^0.5 \ F0;
   W = w;   
   plot(W,A(1),'*');
   grid on;
end