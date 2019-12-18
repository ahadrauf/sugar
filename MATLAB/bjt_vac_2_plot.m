clear all
%%
tic
fileID = fopen('bjt_vac_2_data.txt','r');
formatSpec = '%g %g %g %g';
sizeA = [4 Inf];
A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
M=A';
toc
figure(1); clf; grid on; hold on;
plot( M(:,1)/1000/2/pi, M(:,2)/1e-6,'k' );
plot( M(:,1)/1000/2/pi, M(:,3)/1e-6,'r' );
plot( M(:,1)/1000/2/pi, M(:,4)/1e-6,'b' );
xlabel('Frequency [Hz]')
ylabel('Amplitude [um]')


