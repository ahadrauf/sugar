clear all
%%
tic
fileID = fopen('bjt_k_data.m','r');
formatSpec = '%g %g %g %g';
sizeA = [4 Inf];
A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
M=A';
toc

figure(1); clf; grid on; hold on;
plot( M(:,1)/1000, M(:,2)/1e-6,'r' );
plot( M(:,1)/1000, M(:,3)/1e-6,'b' );
plot( M(:,1)/1000, M(:,4)/1e-6,'g' );
xlabel('Frequency [kHz]')
ylabel('Displacement [um]')
title('g=Ke=-3K/4	b=Ke=0	r=Ke=3K')

%%
tic
fileID = fopen('bjt_m_data.m','r');
formatSpec = '%g %g %g %g';
sizeA = [4 Inf];
A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
M=A';
toc

figure(2); clf; grid on; hold on;
plot( M(:,1)/1000, M(:,2)/1e-6,'r' );
plot( M(:,1)/1000, M(:,3)/1e-6,'b' );
plot( M(:,1)/1000, M(:,4)/1e-6,'g' );
xlabel('Frequency [kHz]')
ylabel('Displacement [um]')
title('g=Me=-3M/4 b=Me=0 r=Me=3M')

%%
tic
fileID = fopen('bjt_d.m','r');
formatSpec = '%g %g %g %g';
sizeA = [4 Inf];
A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
M=A';
toc

figure(3); clf; grid on; hold on;
plot( M(:,1)/1000, M(:,2)/1e-6,'r' );
plot( M(:,1)/1000, M(:,3)/1e-6,'b' );
plot( M(:,1)/1000, M(:,4)/1e-6,'g' );
xlabel('Frequency [kHz]')
ylabel('Displacement [um]')
title('g=De=-D/2 b=De=0 r=De=D')




