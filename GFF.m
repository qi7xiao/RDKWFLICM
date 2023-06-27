function y=GFF(im,sigd,sigr)
% bilateral filter双边滤波器
% 函数输入：
%           im    输入的图像
%           sigd  空间内核的时域参数
%           sigr  内核参数强度变化范围
% 函数输出：
%          out  滤波图像 = output imagespatial kernel
% sigr=(n*100)^2/(.003*(sigd^2));  % 自适应R值，n为高斯噪声强度,n=0.001
[~,~,channel]=size(im);
if channel==1
    y=calGFF(im,sigd,sigr);
elseif channel==3
   y(:,:,1)=calGFF(im(:,:,1),sigd,sigr);
   y(:,:,2)=calGFF(im(:,:,2),sigd,sigr);
   y(:,:,3)=calGFF(im(:,:,3),sigd,sigr);
end

function out=calGFF(im,sigd,sigr)
% 高斯频域滤波器 
        
% Padding 扩展图像的边界，防止滑动窗口边界值溢出
proci=padarray(im,[sigd sigd],'replicate');
[row,clm]=size(proci);    % Size of image
    proci = double(proci)/255;   % 转换为double类型
K=sigd;
L=-K:K;
c=K+1;   % 中心元素位置
for r=(1+K):(row-K)          % 行    
    for s=(1+K):(clm-K)      % 列                     
            win=proci((r+L),(s+L));% 获取窗口   
            I=win; % 灰度矩阵
            win=win(c,c)-win; % 相对中心点处的强度差异，中心点为参考灰度值
            win=sqrt(win.^2); % 保证win中的每一个元素为正
            Gwi=exp(-(win.^2)/(2*(sigr^2))); % 高斯函数      
            weights=(Gwi)/sum(sum(Gwi)); % 高斯权值
            proci(r,s) =sum(sum(weights.*I));               % 得到当前双边滤波值  
    end
end

% 移除边界扩展值
out=im2uint8(proci(K+1:end-K,K+1:end-K));   % 类型转换


