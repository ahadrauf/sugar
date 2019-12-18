clear all
%%
tic
fileID = fopen('bjt_freq_same_Fdr_data.m','r');
formatSpec = '%g %g %g %g';
sizeA = [4 Inf];
A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
M=A';
toc

figure(1); clf; grid on; hold on;
plot( M(:,1)/1e3, M(:,2)/1e-6,'b' );
plot( M(:,1)/1e3, M(:,3)/1e-6,'r' );
plot( M(:,1)/1e3, M(:,4)/1e-6,'g' );
xlabel('Frequency [kHz]')
ylabel('Displacement [um]')

