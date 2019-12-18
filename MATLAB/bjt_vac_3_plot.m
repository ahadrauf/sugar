clear all
%%
tic
fileID = fopen('bjt_vac_3_data.txt','r');
formatSpec = '%g %g %g %g';
sizeA = [4 Inf];
A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
M=A';
toc
figure(1); clf; grid on; hold on;
plot( M(:,1)/1e-3, M(:,2)/1e-6,'b' );
xlabel('time [ms]')
ylabel('Displacement [um]')


