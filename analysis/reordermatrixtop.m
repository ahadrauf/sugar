function [M]=reordermatrixtop(M,node1ids,node2ids)
[M]=reordermatrixdimensions(M,node1ids,[1:6]); %reorder 
[M]=reordermatrixdimensions(M,node2ids,[7:12]); %reorder 

