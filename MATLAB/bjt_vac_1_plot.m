clear all
%%
tic
fileID = fopen('bjt_vac_00D_data.txt','r');
formatSpec = '%g %g %g %g';
sizeA = [4 Inf];
A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
M=A';
toc
figure(1); clf; grid on; hold on;
plot( M(:,1)/1e-3, M(:,2)/1e-6,'k' );
xlabel('Time [ms]')
ylabel('Displacement [um]')

tic
fileID = fopen('bjt_vac_01D_data.txt','r');
formatSpec = '%g %g %g %g';
sizeA = [4 Inf];
A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
M=A';
toc
figure(1); 
plot( M(:,1)/1e-3, M(:,2)/1e-6,'b' );
xlabel('Time [ms]')
ylabel('Displacement [um]')

tic
fileID = fopen('bjt_vac_05D_data.txt','r');
formatSpec = '%g %g %g %g';
sizeA = [4 Inf];
A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
M=A';
toc
figure(1); 
plot( M(:,1)/1e-3, M(:,2)/1e-6,'r' );
xlabel('Time [ms]')
ylabel('Displacement [um]')

tic
fileID = fopen('bjt_vac_D_data.txt','r');
formatSpec = '%g %g %g %g';
sizeA = [4 Inf];
A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
M=A';
toc
figure(1); 
plot( M(:,1)/1e-3, M(:,2)/1e-6,'g' );
xlabel('Time [ms]')
ylabel('Displacement [um]')

figure(1); 
D = 1.55e-7;
gamma = D / 2 / 8e-10;
plot( M(:,1)/1e-3, 10 .* exp( -gamma .* M(:,1)) ,'k' );
D = 1.55e-7;
gamma = 0.1*D / 2 / 8e-10;
plot( M(:,1)/1e-3, 10 .* exp( -gamma .* M(:,1)) ,'k' );
D = 1.55e-7;
gamma = 0.5*D / 2 / 8e-10;
plot( M(:,1)/1e-3, 10 .* exp( -gamma .* M(:,1)) ,'k' );
D = 1.55e-7;
gamma = 0*D / 2 / 8e-10;
plot( M(:,1)/1e-3, 10 .* exp( -gamma .* M(:,1)) ,'k' );
