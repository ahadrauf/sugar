%function [a] = invert(a,nmax,ndm)

a=[             123370.055                         0               -292815.081
                         0                123370.005               -53592.1161
               -292815.081               -53592.1161                2267716.28]
a=eye(3)*3
b =[          0.0743150829             -0.0301989212             -0.0484094359
             -0.0301357926              0.0751907801              0.0311122075
             -0.0483968731               0.031071315               0.071862906];
ia=inv(a)

nmax=3; 
ndm=3; 

    for n = 1 : nmax
        d = a(n,n);        
        for j = 1 : nmax
            a(n,j) = -a(n,j)/d;
        end
        for i = 1 : nmax        
            if (n ~= i)               
                for j = 1 : nmax
                    if(n ~= j) 
                        a(i,j) = a(i,j) + a(i,n)*a(n,j);
                    end
                end
                a(i,n) = a(i,n)/d;
            end
        end
        a(n,n) = 1/d;      
    end %for n      

a1=a

    