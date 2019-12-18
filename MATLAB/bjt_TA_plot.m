clear all
%%
tic
fileID = fopen('bjt_TA_D_data.m','r');
formatSpec = '%g %g %g %g';
sizeA = [4 Inf];
A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
M=A';
toc

figure(1); clf; grid on; hold on;
plot( M(:,1), M(:,2),'b' );
xlabel('Time')
ylabel('Displacement [um]')
title('TA')

tic
fileID = fopen('bjt_TA_Dcrit_data.m','r');
formatSpec = '%g %g %g %g';
sizeA = [4 Inf];
A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
M=A';
toc

plot( M(:,1)+1, M(:,2),'b' );
