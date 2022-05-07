clear all;clc;


x = [448 815 1196 1515 1946 2217 866 1285 1629 2091 2379 866 1285 1629 2091 2379 ]; 
x=x';

%load DWI img
[Filename Filepath] = uigetfile('*.img','select the dwi img');
Datapath = [Filepath filesep Filename];
DWIData = spm_read_vols(spm_vol(Datapath));

%% load mask 
[Filename1 Filepath1] = uigetfile('*.nii','select the mask roi img');
Datapath1 = [Filepath1 filesep Filename1];
Mask = spm_read_vols(spm_vol(Datapath1));
Mask = squeeze(Mask(:,:,1));

%mask threshold dwi data
for bnum = 1:size(DWIData,3)
    Mask_DWIData(:,:,bnum) = DWIData(:,:,bnum).*Mask;
end


%%
bnum_perdir = 5; 

[size1 size2 size3] = size(Mask_DWIData);

%把三个方向的图像分别提取出来
Mask_DWIData_x = zeros(size1,size2,bnum_perdir+1);
Mask_DWIData_y = zeros(size1,size2,bnum_perdir+1);
Mask_DWIData_z = zeros(size1,size2,bnum_perdir+1);

Mask_DWIData_x(:,:,1) = Mask_DWIData(:,:,1);
Mask_DWIData_x(:,:,2:end) = Mask_DWIData(:,:,2:bnum_perdir+1);

Mask_DWIData_y(:,:,1) = Mask_DWIData(:,:,1);
Mask_DWIData_y(:,:,2:end) = Mask_DWIData(:,:,bnum_perdir+2:bnum_perdir*2+1);

Mask_DWIData_z(:,:,1) = Mask_DWIData(:,:,1);
Mask_DWIData_z(:,:,2:end) = Mask_DWIData(:,:,bnum_perdir*2+2:end);

%把三个方向的b值分别提取出来
xx = zeros(bnum_perdir+1,1);
xy = zeros(bnum_perdir+1,1);
xz = zeros(bnum_perdir+1,1);

xx(1,:) = x(1,:);
xx(2:end,:) = x(2:bnum_perdir+1,:);

xy(1,:) = x(1,:);
xy(2:end,:) = x(bnum_perdir+2:bnum_perdir*2+1,:);

xz(1,:) = x(1,:);
xz(2:end,:) = x(bnum_perdir*2+2:end,:);

%% fitting vol by vol 
%第一个方向，x方向
for m=1:size(Mask_DWIData_x,1)
    for n=1:size(Mask_DWIData_x,2)
%         if Mask_DWIData(m,n,1)~=0
        yx = log(squeeze(Mask_DWIData_x(m,n,:)));
        px = polyfit(xx,yx,1);
        ADC_x(m,n) = abs(px(1));
%         fitcurve(m,n) = exp(p(2)+x*p(1));
%         end
    end
end

figure('name','ADC_map_x');
imagesc(ADC_x);colorbar('NorthOutside'); colormap jet; caxis([0.000,0.0015]); %colorbar范围可修改

figure('name','DWI_img_x');
imagesc(Mask_DWIData_x(:,:,size(xx,1))); colormap gray;

%第二个方向，y方向
for m=1:size(Mask_DWIData_y,1)
    for n=1:size(Mask_DWIData_y,2)
%         if Mask_DWIData(m,n,1)~=0
        yy = log(squeeze(Mask_DWIData_y(m,n,:)));
        py = polyfit(xy,yy,1);
        ADC_y(m,n) = abs(py(1));
%         fitcurve(m,n) = exp(p(2)+x*p(1));
%         end
    end
end

figure('name','ADC_map_y');
imagesc(ADC_y);colorbar('NorthOutside'); colormap jet; caxis([0.000,0.0015]); %colorbar范围可修改

figure('name','DWI_img_y');
imagesc(Mask_DWIData_y(:,:,size(xy,1))); colormap gray;

%第三个方向，z方向
for m=1:size(Mask_DWIData_z,1)
    for n=1:size(Mask_DWIData_z,2)
%         if Mask_DWIData(m,n,1)~=0
        yz = log(squeeze(Mask_DWIData_z(m,n,:)));
        pz = polyfit(xz,yz,1);
        ADC_z(m,n) = abs(pz(1));
%         fitcurve(m,n) = exp(p(2)+x*p(1));
%         end
    end
end

figure('name','ADC_map_z');
imagesc(ADC_z);colorbar('NorthOutside'); colormap jet; caxis([0.000,0.0015]); %colorbar范围可修改

figure('name','DWI_img_z');
imagesc(Mask_DWIData_z(:,:,size(xz,1))); colormap gray;

ADC_XYZ=(ADC_x+ADC_y+ADC_z)/3;
figure('name','ADC_map');
imagesc(ADC_XYZ);colorbar('WestOutside'); colormap jet; caxis([0.000,0.0015]); %colorbar范围可修改

%% ROIfitting_curve

figure('name','b0_img');
imagesc(DWIData(:,:,1)); colormap gray;
hold on;
ROISel = roipoly;

SelctSig_x = [];
    for i = 1:size(Mask_DWIData_x,3)
        ROValue_x    = Mask_DWIData_x(:,:,i).*ROISel;
        ROILoc_x     = find(ROValue_x~=0);
        ROIValue_x   = ROValue_x(ROILoc_x);
        SelctSig_x   = [SelctSig_x, ROIValue_x];
    end   
    ROImeanValue_x  = mean(SelctSig_x);
    ROImeanValue_x = (squeeze(ROImeanValue_x))';
    
 SelctSig_y = [];
    for i = 1:size(Mask_DWIData_y,3)
        ROValue_y    = Mask_DWIData_y(:,:,i).*ROISel;
        ROILoc_y     = find(ROValue_y~=0);
        ROIValue_y   = ROValue_y(ROILoc_y);
        SelctSig_y   = [SelctSig_y, ROIValue_y];
    end   
    ROImeanValue_y  = mean(SelctSig_y);
    ROImeanValue_y = (squeeze(ROImeanValue_y))'; 
    
 SelctSig_z = [];
    for i = 1:size(Mask_DWIData_z,3)
        ROValue_z    = Mask_DWIData_z(:,:,i).*ROISel;
        ROILoc_z     = find(ROValue_z~=0);
        ROIValue_z   = ROValue_z(ROILoc_z);
        SelctSig_z   = [SelctSig_z, ROIValue_z];
    end   
    ROImeanValue_z  = mean(SelctSig_z);
    ROImeanValue_z = (squeeze(ROImeanValue_z))'; 
    
    %% 画圆形ROI

figure('name','b0_img');
imagesc(DWIData(:,:,1)); colormap gray;
hold on;
% ROISel = roipoly;
    
[imgW,imgH] = size(squeeze(DWIData(:,:,2)));
t = linspace(0, 2*pi, 50);   %# approximate circle with 50 points
r = 2;                      %半径
[x0,y0]= ginput(1);
c = [x0 y0];
% c = [157 107];               %圆心坐标
%get circular mask
ROISel = poly2mask(r*cos(t)+c(1), r*sin(t)+c(2), imgW, imgH);
hold on;
plot(r*cos(t)+c(1), r*sin(t)+c(2));

SelctSig_x = [];
    for i = 1:size(Mask_DWIData_x,3)
        ROValue_x    = Mask_DWIData_x(:,:,i).*ROISel;
        ROILoc_x     = find(ROValue_x~=0);
        ROIValue_x   = ROValue_x(ROILoc_x);
        SelctSig_x   = [SelctSig_x, ROIValue_x];
    end   
    ROImeanValue_x  = mean(SelctSig_x);
    ROImeanValue_x = (squeeze(ROImeanValue_x))';
    
 SelctSig_y = [];
    for i = 1:size(Mask_DWIData_y,3)
        ROValue_y    = Mask_DWIData_y(:,:,i).*ROISel;
        ROILoc_y     = find(ROValue_y~=0);
        ROIValue_y   = ROValue_y(ROILoc_y);
        SelctSig_y   = [SelctSig_y, ROIValue_y];
    end   
    ROImeanValue_y  = mean(SelctSig_y);
    ROImeanValue_y = (squeeze(ROImeanValue_y))'; 
    
 SelctSig_z = [];
    for i = 1:size(Mask_DWIData_z,3)
        ROValue_z    = Mask_DWIData_z(:,:,i).*ROISel;
        ROILoc_z     = find(ROValue_z~=0);
        ROIValue_z   = ROValue_z(ROILoc_z);
        SelctSig_z   = [SelctSig_z, ROIValue_z];
    end   
    ROImeanValue_z  = mean(SelctSig_z);
    ROImeanValue_z = (squeeze(ROImeanValue_z))'; 
    
    
    %% 后续可用cftool拟合，得到方程后，画出曲线和点图
    x1 = linspace(0,2500,25000); 
    expfitcurve = 90.56*exp(-0.0005698*x1); %f(x)=a*exp(b*x); -b为ADC值，需要从cftool里获得
    figure('name','ROI_plot')
    plot(xx,ROImeanValue_x,'ro',x1,expfitcurve,'-');
    
   %% 线性拟合
    y1 = log(ROImeanValue);
    [p1 S]= polyfit(x,y1,1); %S中得不到残差？
   
    x1 = linspace(0,2500,25000); %范围及点数是手动设置的
    logy = polyval(p1,x1);
    fitcurve = exp(logy);
    
%     A = exp(p1(2));% y=A*exp(-b*D) %-p1(1)为ADC
%     fitcurve = A*exp(x1*p1(1));
    
    figure('name','ROI_plot')
    plot(x,ROImeanValue,'ro',x1,fitcurve,'-');
        
    %% 指数拟合
    myfun = @(xexp,x) xexp(1)*exp(xexp(2)*x);
    xexp0= [87.5561714187041 -0.000541073595658469]; %starting guess  %cftool代码中有，cftool界面 fit options按钮中有
    [xexp,resnorm] = lsqcurvefit(myfun,xexp0,x,ROImeanValue);
    
    x1 = linspace(0,2500,25000); 
    expfitcurve = xexp(1)*exp(xexp(2)*x1); %-xexp(2)为ADC值
    figure('name','ROI_plot')
    plot(x,ROImeanValue,'ro',x1,expfitcurve,'-');
    
%% cftool代码
    ft = fittype( 'exp1' ); %f(x)=a*exp(b*x)
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Algorithm = 'Levenberg-Marquardt';
    opts.Display = 'Off';
    opts.StartPoint = [87.5561714187041 -0.000541073595658469]; %cftool可以自动给

    % Fit model to data.
    [fitresult, gof,output] = fit( x, ROImeanValue, ft, opts );
    gof
    R2 = gof.rsquare;
    SSE = gof.sse;
    
    % Plot fit with data.
    figure( 'Name', 'exp fit' );
    h = plot( fitresult, x, ROImeanValue );
    legend( h, 'ROImeanValue vs. x', 'exp fit', 'Location', 'NorthEast' );
    % Label axes
    xlabel( 'x' );
    ylabel( 'ROImeanValue' );



