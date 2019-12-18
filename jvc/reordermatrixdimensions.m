%This function reorders the index of an N-dimensional N*N matrix M or N*1 vector.
%Eg: reorder={[1 2 3],[7 8 9]} will switch [1 2 3] components to [7 8 9] and vise versa.
function [N]=reordermatrixdimensions(M,a,b)
%By Jason Vaughn Clark - November2002
N=M;
%switch rows
N(a,:)=M(b,:);
N(b,:)=M(a,:);
M=N;
%switch columns
if min(size(M))>1
   N(:,a)=M(:,b);
   N(:,b)=M(:,a);
end

