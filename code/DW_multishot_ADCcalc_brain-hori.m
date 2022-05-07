clear all;clc;

% load DWI img（4d horizontal）
[Filename Filepath] = uigetfile('*.nii','select the dwi img');
Datapath = [Filepath filesep Filename];
DWIData = spm_read_vols(spm_vol(Datapath));

[Filename1 Filepath1] = uigetfile('*.nii','select the mask img');
Datapath1 = [Filepath1 filesep Filename1];
MaskData = spm_read_vols(spm_vol(Datapath1));

% Slicenum =10;

b_value = [168 1575 1820 1820 ];


%提取b=0
DWIData_b0 = DWIData(:,:,:,1) ;
% DWIData_b0_brain = DWIData(:,:,:,1).*MaskData ;


%% 提取出x方向
DWIData_b1000_x = DWIData(:,:,:,2) ;

%计算ADC值
[size_x,size_y,size_z]=size(DWIData_b1000_x);
for nslice =1:size_z
    for i=1:size_x
      for j=1:size_y
      ADC_x(i,j,nslice)= (log(DWIData_b1000_x(i,j,nslice)/DWIData_b0(i,j,nslice)))/(b_value(1)-b_value(2));
      end
    end
end

ADC_x_brain = ADC_x.*MaskData;

figure('name','ADC_x_img');imagesc(ADC_x_brain(:,:,4)',[0.0000,0.0013]);colormap(hot(30));
figurenamex = [Filepath,'ADC_x'];
print(gcf,'-dtiff' ,'-r300',figurenamex);

%% 提取出y方向
DWIData_b1000_y = DWIData(:,:,:,3) ;

%计算ADC值
for nslice =1:size_z
    for i=1:size_x
      for j=1:size_y
      ADC_y(i,j,nslice)= (log(DWIData_b1000_y(i,j,nslice)/DWIData_b0(i,j,nslice)))/(b_value(1)-b_value(3));
      end
    end
end

ADC_y_brain = ADC_y.*MaskData;

figure('name','ADC_y_img');imagesc(ADC_y_brain(:,:,4)',[0.0000,0.0013]);colormap(hot(30));
figurenamey = [Filepath,'ADC_y'];
print(gcf,'-dtiff' ,'-r300',figurenamey);


%% 提取出z方向
DWIData_b1000_z = DWIData(:,:,:,4) ;

%计算ADC值
for nslice =1:size_z
    for i=1:size_x
      for j=1:size_y
      ADC_z(i,j,nslice)= (log(DWIData_b1000_z(i,j,nslice)/DWIData_b0(i,j,nslice)))/(b_value(1)-b_value(4));
      end
    end
end

ADC_z_brain = ADC_z.*MaskData;

figure('name','ADC_z_img');imagesc(ADC_z_brain(:,:,4)',[0.0000,0.0013]);colormap(hot(30));
figurenamez = [Filepath,'ADC_z'];
print(gcf,'-dtiff' ,'-r300',figurenamez);


%%  三个方向ADC值进行平均
% ADC=(ADC_x+ADC_y+ADC_z)/3;
ADC=(ADC_x_brain+ADC_y_brain+ADC_z_brain)/3;

figure('name','ADC_img');imagesc(ADC(:,:,4)',[0.0000,0.0013]);colormap(hot(30));
figurename = [Filepath,'ADC'];
print(gcf,'-dtiff' ,'-r300',figurename);

%% 读一个3d的头文件
[Filename2 Filepath2] = uigetfile('*.nii','select the header img');
Datapath2 = [Filepath2 filesep Filename2];
Header = spm_vol(Datapath2);


%% 保存ADC为img图像
[FileName3, FileDir3]         = uiputfile('*.img','ADC name');
FilePath3                    = [FileDir3 FileName3];
RestHead1                  = Header;
% RestHead1.dim(1)           = 128;
% RestHead1.dim(2)           = 100;
% RestHead1.mat(1,1:3)=[-0.140625000000000,0,0];
% RestHead1.mat(2,1:3)=[0,0.18,0];
% RestHead1.mat(3,1:3)=[0,0,0.800000011920929];
% RestHead1.mat(1:4,4)=[9.07031250000000;-9.10000013560057;-4.8;1];
RestHead1.dt=[16,0];
RestHead1.fname             = FilePath3;
RestHead1                   = spm_create_vol(RestHead1);
spm_write_plane(RestHead1,ADC,':');

%% 保存ADC_x为img图像
[FileName3, FileDir3]         = uiputfile('*.img','ADC name');
FilePath3                    = [FileDir3 FileName3];
RestHead1                  = Header;
% RestHead1.dim(1)           = 128;
% RestHead1.dim(2)           = 100;
% RestHead1.mat(1,1:3)=[-0.140625000000000,0,0];
% RestHead1.mat(2,1:3)=[0,0.18,0];
% RestHead1.mat(3,1:3)=[0,0,0.800000011920929];
% RestHead1.mat(1:4,4)=[9.07031250000000;-9.10000013560057;-4.8;1];
RestHead1.dt=[16,0];
RestHead1.fname             = FilePath3;
RestHead1                   = spm_create_vol(RestHead1);
spm_write_plane(RestHead1,ADC_x,':');


%% 保存ADC_y为img图像
[FileName3, FileDir3]         = uiputfile('*.img','ADC name');
FilePath3                    = [FileDir3 FileName3];
RestHead1                  = Header;
% RestHead1.dim(1)           = 128;
% RestHead1.dim(2)           = 100;
% RestHead1.mat(1,1:3)=[-0.140625000000000,0,0];
% RestHead1.mat(2,1:3)=[0,0.18,0];
% RestHead1.mat(3,1:3)=[0,0,0.800000011920929];
% RestHead1.mat(1:4,4)=[9.07031250000000;-9.10000013560057;-4.8;1];
RestHead1.dt=[16,0];
RestHead1.fname             = FilePath3;
RestHead1                   = spm_create_vol(RestHead1);
spm_write_plane(RestHead1,ADC_y_brain,':');


%% 保存ADC_z为img图像
[FileName3, FileDir3]         = uiputfile('*.img','ADC name');
FilePath3                    = [FileDir3 FileName3];
RestHead1                  = Header;
% RestHead1.dim(1)           = 128;
% RestHead1.dim(2)           = 100;
% RestHead1.mat(1,1:3)=[-0.140625000000000,0,0];
% RestHead1.mat(2,1:3)=[0,0.18,0];
% RestHead1.mat(3,1:3)=[0,0,0.800000011920929];
% RestHead1.mat(1:4,4)=[9.07031250000000;-9.10000013560057;-4.8;1];
RestHead1.dt=[16,0];
RestHead1.fname             = FilePath3;
RestHead1                   = spm_create_vol(RestHead1);
spm_write_plane(RestHead1,ADC_z,':');


%% 保存b0图为img
[FileName4, FileDir4]         = uiputfile('*.img','b0 name');
FilePath4                    = [FileDir4 FileName4];
RestHead2                  = Header;
% RestHead2.dim(1)           = 128;
% RestHead2.dim(2)           = 100;
% RestHead2.mat(1,1:3)=[-0.140625000000000,0,0];
% RestHead2.mat(2,1:3)=[0,0.18,0];
% RestHead2.mat(3,1:3)=[0,0,0.800000011920929];
% RestHead2.mat(1:4,4)=[9.07031250000000;-9.10000013560057;-4.8;1];
RestHead1.dt=[16,0];
RestHead2.fname             = FilePath4;
RestHead2                  = spm_create_vol(RestHead2);
spm_write_plane(RestHead2,DWIData_b0_brain,':');








