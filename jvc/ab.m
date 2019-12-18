function ab

fid=fopen('cif.txt', 'w')    

a = 3;
b = 41;

fprintf(fid,'B %d %d \n',a,b);
fclose(fid)  