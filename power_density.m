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
TOT_POWER=(160*10^6)/37;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DEFINIÇAO DO MALHA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx=length(DET1X);
Ny=length(DET1Y);
Nz=length(DET1Z);

DX=DET1X(end,1)-DET1X(end-1,1);
DY=DET1Y(end,1)-DET1Y(end-1,1);
DZ=DET1Z(end,1)-DET1Z(end-1,1);

for i=1:length(DET1X)
    dx(i)=DET1X(i+1)-DET1X(i);
    dy(i)=DET1Y(i+1)-DET1Y(i);
end

for i=1:length(DET1Z)
    dz(i)=DET1Z(i+1)-DET1Z(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VETORES AXIAIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xa=[1:Nx]*DX;
Ya=[1:Ny]*DY;
Za=DET1Z(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MALHA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y,Z]=meshgrid(Xa,Ya,Za);

AX=sum(DET1X);
AY=sum(DET1Y);
AZ=sum(DET1Z);

LX=-AX(1)+AX(2);
LY=-AY(1)+AY(2);
LZ=-AZ(1)+AZ(2);

for i=1:length(DET1X)
    Xx(i)=DET1X(i+1)-DET1X(i);
    Yy(i)=DET1Y(i+1)-DET1Y(i);
end

for i=1:length(DET1Z)
    Zz(i)=DET1Z(i+1)-DET1Z(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CALCULO DO VOLUME DA MALHA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(DET1X)
    for j=1:length(DET1Y)
        for k=1:length(DET1Z)
        VOLUME_MESH(i,j,k)=Xa(i)*Ya(j)*Za(k);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NORMALIZAÇÃO DOS RESULTADOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TOT=sum(DET1(:,11));
REACT=max(DET1(:,11));
REAC=(DET1(1:end,11));

REAC_X=DET1(1:length(DET1Y)*length(DET1Z):end,11);
REAC_Y=DET1(11:length(DET1X)*length(DET1Z):end,11);
REAC_Z=DET1(length(DET1X)*length(DET1Y):length(DET1X)*length(DET1Y):end,11);

TOT_REACTION=sum(DET1(:,11))

PDENS_Z=(((DET1(:,11))/TOT_REACTION)*TOT_POWER)/(LZ/Nz);
PDENS_X=(((DET1(:,11))/TOT_REACTION)*TOT_POWER)/(LX*LZ/Nx*Nz);
PDENS_Y=(((DET1(:,11))/TOT_REACTION)*TOT_POWER)/(LY*LZ/Ny*Nz);

fprintf('Power Density Y-axis[W/cm] %d \n',max(PDENS_Y));
fprintf('Power Density X-axis[W/cm] %d \n',max(PDENS_X));
fprintf('Power Density Z-axis[W/cm] %d \n',max(PDENS_Z));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CONVERSÃO DE UMA CÉLULA DE DENSIDADE DE POTÊNCIA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PDENSZZ=cell(1,length(DET1X)*length(DET1Y));
for i=1:length(DET1X)*length(DET1Y);
    PDENSZZ{i}=PDENS_Z(i:length(DET1X)*length(DET1Y):end,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CONVERSÃO DE UMA CÉLULA DE MATRIZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAT_PWR=cell2mat(PDENSZZ);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MÁX DA DENSIDADE DE POTÊNCIA EM CADA LINHA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(DET1Z);
    MP(i,:)=max(MAT_PWR(i,:));
end
    hold on
    figure (1)
    hold on
    plot(Za,MP,'.')
    title('Figura - 1')
    xlabel('Z axial[cm]')
    ylabel('Densidade de Potência[W/cm]')
    
    hold on
    figure (2)
    hold on
    bar(MP)
    title('Figura - 2')
    xlabel('Z Grid')
    ylabel('Densidade de Potência[W/cm]')
    
qz=[MP(1:1:end)];
zz=[DET1Z(:,1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GERAR UM ARQUIVO TXT COM RESULTADOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID1 = fopen('z_axis_nuscale.txt','w');
fprintf(fileID1,'%6s \n','Z[m]');
fprintf(fileID1,'%f  \n',zz);
fileID2 = fopen('qeoc_nuscale.txt','w');
fprintf(fileID2,'%6s \n','Q[W/cm]');
fprintf(fileID2,'%f  \n',qz);
fclose(fileID1);
fclose(fileID2);


