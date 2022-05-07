
clear all;clc;

% load DWI img
[Filename Filepath] = uigetfile('*.img','select the dwi img');
Datapath = [Filepath filesep Filename];
DWIData = spm_read_vols(spm_vol(Datapath));

[Filename1 Filepath1] = uigetfile('*.nii','select the mask img');
Datapath1 = [Filepath1 filesep Filename1];
MaskData = spm_read_vols(spm_vol(Datapath1));

Slicenum =12;

b_value = [168 1466 1699 1699 ];


DWIData_b0 = DWIData(:,:,1:Slicenum) ;


%% 提取出x方向
DWIData_b1000_x = DWIData(:,:,1+Slicenum:Slicenum*2) ;

%计算ADC值
[size_x,size_y,size_z]=size(DWIData_b1000_x);
for nslice =1:size_z
    for i=1:size_x
      for j=1:size_y
      ADC_x(i,j,nslice)= (log(DWIData_b1000_x(i,j,nslice)/DWIData_b0(i,j,nslice)))/(b_value(1)-b_value(2));
      end
    end
end

ADC_x_brain = ADC_x.*MaskData(:,:,1:Slicenum);

figure('name','ADC_x_img');imagesc(ADC_x_brain(:,:,9)',[0.0000,0.0011]);colormap(hot(30));
figurenamex = [Filepath,'ADC_x'];
print(gcf,'-dtiff' ,'-r300',figurenamex);


%% 提取出y方向
DWIData_b1000_y = DWIData(:,:,1+Slicenum*2:Slicenum*3) ;

%计算ADC值
for nslice =1:size_z
    for i=1:size_x
      for j=1:size_y
      ADC_y(i,j,nslice)= (log(DWIData_b1000_y(i,j,nslice)/DWIData_b0(i,j,nslice)))/(b_value(1)-b_value(3));
      end
    end
end

ADC_y_brain = ADC_y.*MaskData(:,:,1:Slicenum);

figure('name','ADC_y_img');imagesc(ADC_y_brain(:,:,9)',[0.0000,0.0011]);colormap(hot(30));
figurenamex = [Filepath,'ADC_y'];
print(gcf,'-dtiff' ,'-r300',figurenamex);

%% 提取出z方向
DWIData_b1000_z = DWIData(:,:,1+Slicenum*3:Slicenum*4) ;

%计算ADC值
for nslice =1:size_z
    for i=1:size_x
      for j=1:size_y
      ADC_z(i,j,nslice)= (log(DWIData_b1000_z(i,j,nslice)/DWIData_b0(i,j,nslice)))/(b_value(1)-b_value(4));
      end
    end
end

ADC_z_brain = ADC_z.*MaskData(:,:,1:Slicenum);

figure('name','ADC_z_img');imagesc(ADC_z_brain(:,:,9)',[0.0000,0.0011]);colormap(hot(30));
figurenamex = [Filepath,'ADC_z'];
print(gcf,'-dtiff' ,'-r300',figurenamex);


%% 三个方向ADC值进行平均
ADC=(ADC_x_brain+ADC_y_brain+ADC_z_brain)/3;

figure('name','ADC_img');imagesc(ADC(:,:,9)',[0.0000,0.0011]);colormap(hot(30));
figurenamex = [Filepath,'ADC'];
print(gcf,'-dtiff' ,'-r300',figurenamex);


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
RestHead1.dim(1)           = 128;
RestHead1.dim(2)           = 100;
RestHead1.mat(1,1:3)=[-0.140625000000000,0,0];
RestHead1.mat(2,1:3)=[0,0.18,0];
RestHead1.mat(3,1:3)=[0,0,0.800000011920929];
RestHead1.mat(1:4,4)=[9.07031250000000;-9.10000013560057;-4.8;1];
RestHead1.dt=[16,0];
RestHead1.fname             = FilePath3;
RestHead1                   = spm_create_vol(RestHead1);
spm_write_plane(RestHead1,ADC_x_brain,':');


%% 保存ADC_y为img图像
[FileName3, FileDir3]         = uiputfile('*.img','ADC name');
FilePath3                    = [FileDir3 FileName3];
RestHead1                  = Header;
RestHead1.dim(1)           = 128;
RestHead1.dim(2)           = 100;
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
RestHead1.dim(1)           = 128;
RestHead1.dim(2)           = 100;
RestHead1.mat(1,1:3)=[-0.140625000000000,0,0];
RestHead1.mat(2,1:3)=[0,0.18,0];
RestHead1.mat(3,1:3)=[0,0,0.800000011920929];
RestHead1.mat(1:4,4)=[9.07031250000000;-9.10000013560057;-4.8;1];
RestHead1.dt=[16,0];
RestHead1.fname             = FilePath3;
RestHead1                   = spm_create_vol(RestHead1);
spm_write_plane(RestHead1,ADC_z_brain,':');




%% 保存b0图为img
[FileName4, FileDir4]         = uiputfile('*.img','b0 name');
FilePath4                    = [FileDir4 FileName4];
RestHead2                  = Header;
RestHead2.dim(1)           = 128;
RestHead2.dim(2)           = 100;
% RestHead2.mat(1,1:3)=[-0.140625000000000,0,0];
% RestHead2.mat(2,1:3)=[0,0.18,0];
% RestHead2.mat(3,1:3)=[0,0,0.800000011920929];
% RestHead2.mat(1:4,4)=[9.07031250000000;-9.10000013560057;-4.8;1];
RestHead1.dt=[4,0];
RestHead2.fname             = FilePath4;
RestHead2                  = spm_create_vol(RestHead2);
DWIData_b0_brain = DWIData_b0.*MaskData(:,:,1:Slicenum);
spm_write_plane(RestHead2,DWIData_b0_brain,':');








