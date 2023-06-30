
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program Power Density
% Author : Diego Manoel 
%
% Insert [file.det] and Run
% The file have DET1 ; DET1X;DET1Y;DET1Z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input SERPENT.VTT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Increase counter:
TOT_POWER=(160/37)*10^6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFINIÇAO DO MALHA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nx=length(DET1X);
Ny=length(DET1Y);
Nz=length(DET1Z);

Nx2=length(DET2X);
Ny2=length(DET2Y);
Nz2=length(DET2Z);

Nx3=length(DET3X);
Ny3=length(DET3Y);
Nz3=length(DET3Z);

Nx4=length(DET4X);
Ny4=length(DET4Y);
Nz4=length(DET4Z);

Nx5=length(DET5X);
Ny5=length(DET5Y);
Nz5=length(DET5Z);

DX=DET1X(end,1)-DET1X(end-1,1);
DY=DET1Y(end,1)-DET1Y(end-1,1);
DZ=DET1Z(end,1)-DET1Z(end-1,1);

DX2=DET2X(end,1)-DET2X(end-1,1);
DY2=DET2Y(end,1)-DET2Y(end-1,1);
DZ2=DET2Z(end,1)-DET2Z(end-1,1);


DX3=DET3X(end,1)-DET3X(end-1,1);
DY3=DET3Y(end,1)-DET3Y(end-1,1);
DZ3=DET3Z(end,1)-DET3Z(end-1,1);

DX4=DET4X(end,1)-DET4X(end-1,1);
DY4=DET4Y(end,1)-DET4Y(end-1,1);
DZ4=DET4Z(end,1)-DET4Z(end-1,1);

DX5=DET5X(end,1)-DET5X(end-1,1);
DY5=DET5Y(end,1)-DET5Y(end-1,1);
DZ5=DET5Z(end,1)-DET5Z(end-1,1);
for i=1:length(DET1X)
    dx(i)=DET1X(i+1)-DET1X(i);
    dy(i)=DET1Y(i+1)-DET1Y(i);
    
    dx2(i)=DET2X(i+1)-DET2X(i);
    dy2(i)=DET2Y(i+1)-DET2Y(i);
    
    dx3(i)=DET3X(i+1)-DET3X(i);
    dy3(i)=DET3Y(i+1)-DET3Y(i);
    
    dx4(i)=DET4X(i+1)-DET4X(i);
    dy4(i)=DET4Y(i+1)-DET4Y(i);
    
    dx5(i)=DET5X(i+1)-DET5X(i);
    dy5(i)=DET5Y(i+1)-DET5Y(i);
end

for i=1:length(DET1Z)
    dz(i)=DET1Z(i+1)-DET1Z(i);
    dz2(i)=DET2Z(i+1)-DET2Z(i);
    dz3(i)=DET3Z(i+1)-DET3Z(i);
    dz4(i)=DET4Z(i+1)-DET4Z(i);
    dz5(i)=DET5Z(i+1)-DET5Z(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VETORES AXIAIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xa=[1:Nx]*DX;
Ya=[1:Ny]*DY;
Za=[1:Nz]*DZ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MALHA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y,Z]=meshgrid(Xa,Ya,Za);

AX=sum(DET1X);
AY=sum(DET1Y);
AZ=sum(DET1Z);

AX2=sum(DET2X);
AY2=sum(DET2Y);
AZ2=sum(DET2Z);

AX3=sum(DET3X);
AY3=sum(DET3Y);
AZ3=sum(DET3Z);

AX4=sum(DET4X);
AY4=sum(DET4Y);
AZ4=sum(DET4Z);

AX5=sum(DET5X);
AY5=sum(DET5Y);
AZ5=sum(DET5Z);

LX=-AX(1)+AX(2);
LY=-AY(1)+AY(2);
LZ=-AZ(1)+AZ(2);

LX2=-AX2(1)+AX2(2);
LY2=-AY2(1)+AY2(2);
LZ2=-AZ2(1)+AZ2(2);

LX3=-AX3(1)+AX3(2);
LY3=-AY3(1)+AY3(2);
LZ3=-AZ3(1)+AZ3(2);

LX4=-AX4(1)+AX4(2);
LY4=-AY4(1)+AY4(2);
LZ4=-AZ4(1)+AZ4(2);

LX5=-AX5(1)+AX5(2);
LY5=-AY5(1)+AY5(2);
LZ5=-AZ5(1)+AZ5(2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NORMALIZAÇÃO DOS RESULTADOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TOT=max(DET1(:,11));

REAC=(DET1(1:end,11));

REAC_X=DET1(1:length(DET1Y)*length(DET1Z):end,11);
REAC_Y=DET1(11:length(DET1X)*length(DET1Z):end,11);
REAC_Z=DET1(length(DET1X)*length(DET1Y):length(DET1X)*length(DET1Y):end,11);

TOT_REACTION=sum(DET1(:,11))+sum(DET2(:,11))+sum(DET3(:,11))+sum(DET4(:,11))+sum(DET5(:,11));

PDENS_Z=DET1(:,11)*(TOT_POWER/TOT_REACTION)/(LZ/Nz);

PDENS_Z2=DET2(:,11)*(TOT_POWER/TOT_REACTION)/(LZ2/Nz2);

PDENS_Z3=DET3(:,11)*(TOT_POWER/TOT_REACTION)/(LZ3/Nz3);

PDENS_Z4=DET4(:,11)*(TOT_POWER/TOT_REACTION)/(LZ4/Nz4);

PDENS_Z5=DET5(:,11)*(TOT_POWER/TOT_REACTION)/(LZ5/Nz5);

PDENS_X=(((DET1(:,11))/TOT_REACTION)*TOT_POWER)/(LX*LZ/Nx*Nz);
PDENS_Y=(((DET1(:,11))/TOT_REACTION)*TOT_POWER)/(LY*LZ/Ny*Nz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CONVERSÃO DE UMA CÉLULA DE DENSIDADE DE POTÊNCIA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PDENSZZ=cell(1,length(DET1X)*length(DET1Y));
PDENSZZ2=cell(1,length(DET2X)*length(DET2Y));
PDENSZZ3=cell(1,length(DET3X)*length(DET3Y));
PDENSZZ4=cell(1,length(DET4X)*length(DET4Y));
PDENSZZ5=cell(1,length(DET5X)*length(DET5Y));

for i=1:length(DET1X)*length(DET1Y);
    PDENSZZ{i}=PDENS_Z(i:length(DET1X)*length(DET1Y):end,1);
    PDENSZZ2{i}=PDENS_Z2(i:length(DET2X)*length(DET2Y):end,1);
    PDENSZZ3{i}=PDENS_Z3(i:length(DET3X)*length(DET3Y):end,1);
    PDENSZZ4{i}=PDENS_Z4(i:length(DET4X)*length(DET4Y):end,1);
    PDENSZZ5{i}=PDENS_Z5(i:length(DET5X)*length(DET5Y):end,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CONVERSÃO DE UMA CÉLULA DE MATRIZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAT_PWR=cell2mat(PDENSZZ);
MAT_PWR2=cell2mat(PDENSZZ2);
MAT_PWR3=cell2mat(PDENSZZ3);
MAT_PWR4=cell2mat(PDENSZZ4);
MAT_PWR5=cell2mat(PDENSZZ5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MÁX DA DENSIDADE DE POTÊNCIA EM CADA LINHA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(DET1Z);
    MP(i,:)=max(MAT_PWR(i,:));
    MP2(i,:)=max(MAT_PWR2(i,:));
    MP3(i,:)=max(MAT_PWR3(i,:));
    MP4(i,:)=max(MAT_PWR4(i,:));
    MP5(i,:)=max(MAT_PWR5(i,:));
end
    figure (1)
    hold on
    plot(DET5Z(:,1),MP5,'.')
    title('Seed')
    xlabel('Z axial[cm]')
    ylabel('Densidade de Potência[W/cm]')
    
    figure (2)
    hold on
    plot(DET2Z(:,1),MP2,'.')
    title('Blanket')
    xlabel('Z axial[cm]')
    ylabel('Densidade de Potência[W/cm]')
    
    figure (3)
    hold on
    bar(MP5)
    title('Figura-2')
    xlabel('Z Grid')
    ylabel('Densidade de Potência[W/cm]')

 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GERAR UM ARQUIVO TXT COM RESULTADOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qz1=[MP5(1:1:end)];
zz1=[DET5Z(:,1)];
qz2=[MP2(1:1:end)];
zz2=[DET2Z(:,1)];
qz3=[MP2(1:1:end)+MP(1:1:end)];
zz3=[DET3Z(:,1)];

fileID1 = fopen('zaxis.txt','w');
fprintf(fileID1,'%6s \n','Z[m]');
fprintf(fileID1,'%f  \n',zz1);
fileID2 = fopen('qeoc_sbu.txt','w');
fprintf(fileID2,'%6s \n','Q[W/cm]');
fprintf(fileID2,'%f  \n',qz1);
fileID3 = fopen('zaxis.txt','w');
fprintf(fileID3,'%6s \n','Z[m]');
fprintf(fileID3,'%f  \n',zz2);
fileID3 = fopen('qeoc_sbu2.txt','w');
fprintf(fileID3,'%6s \n','Q[W/cm]');
fprintf(fileID3,'%f  \n',qz2);
fileID4 = fopen('zaxis.txt','w');
fprintf(fileID4,'%6s \n','Z[m]');
fprintf(fileID4,'%f  \n',zz3);
fileID5 = fopen('qeoc_sbu3.txt','w');
fprintf(fileID5,'%6s \n','Q[W/cm]');
fprintf(fileID5,'%f  \n',qz3);
fclose(fileID1);
fclose(fileID2);
fclose(fileID3);
fclose(fileID4);
fclose(fileID5);
