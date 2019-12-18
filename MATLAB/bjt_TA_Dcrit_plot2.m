clear all
%%
tic
fileID = fopen('bjt_TA_Dcrit_data2.m','r');
formatSpec = '%g %g %g %g';
sizeA = [4 Inf];
A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
M=A';
toc

figure(1); clf; grid on; hold on;
plot( M(:,1), 2*M(:,2)/1e-10,'b' );
xlabel('Time [ms]')
ylabel('Displacement [angstrom]')

