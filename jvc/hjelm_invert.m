      function [a,aa]=invert(a,nmax,ndm)
%C   *---------------------------------------------------------------*
%C   |     Invert matrix A (nmax,nmax) when array dimension is NDM   |
%C   *---------------------------------------------------------------*
%      implicit double precision (a-h,o-z)
%      dimension a(ndm,ndm)

      for n = 1:nmax %do 200 n = 1,nmax
        d = a(n*(ndm-1)+n); %d = a(n,n);
        for j = 1:nmax   %do 100 j = 1,nmax
          a(n*(ndm-1)+j) = -a(n*(ndm-1)+j)/d;   %a(n,j) = -a(n,j)/d;
        end   %100   continue
        for i = 1:nmax   %do 150 i = 1,nmax
          if(n~=i) %if(n.eq.i) go to 150
             for j = 1:nmax %do 140 j = 1,nmax
               if(n~=j) a(i*(ndm-1)+j) = a(i*(ndm-1)+j) + a(i*(ndm-1)+n)*a(n*(ndm-1)+j); end %if(n.ne.j) a(i,j) = a(i,j) + a(i,n)*a(n,j)
             end  %140     continue
             a(i*(ndm-1)+n) = a(i*(ndm-1)+n)/d;
           end
         end   %150   continue
         a(n*(ndm-1)+n) = 1/d;  %a(n,n) = 1.d0/d
       end   %200 continue
%      return
%      end

aa(1,1:3)=a(1:3);
aa(2,1:3)=a(4:6);
aa(3,1:3)=a(7:9);
