function [dm,dq,dp] = applied(z,clf,n,dmo,dqo,dpo)  %Evaluate the distributed load functions at point z	 
    %nominal values of the applied forces at point z
    f = 1 - z; 
    switch n            
        case 1  %Compute total transverse loads
            dm = dmo(1) + (dmo(2) + dmo(3)*f)*clf;    
            dq = dqo(1) + (dqo(2) + dqo(3)*f)*clf;
            dp = dpo(1) + (dpo(2) + dpo(3)*f)*clf;        
        case 2  %transverse loads associated with load weight_factor only
            dm = dmo(2) + dmo(3)*f;   
            dq = dqo(2) + dqo(3)*f;
            dp = dpo(2) + dpo(3)*f;
    end
