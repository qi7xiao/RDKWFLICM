clc,clear,close all  % 清理命令区、清理工作区、关闭显示图形
 warning off       % 消除警告
 feature jit off      % 加速代码运行
 
G=imread('D:\总图库\小目标图\ct31.jpg');
tp2=imread('result/tp2.png');

[~,~,rr]=size(G);
if rr>1
    G=rgb2gray(G);
end
im=double(G);

h=50;                                     % 高斯噪声比准差
NI = randn(size(im))*h;                 % 白色高斯噪声
I=im+ abs(NI);   % 把噪声加到原图上面
% I=imnoise(G,'gaussian',0,0.1);
% I=imnoise(G,'salt & pepper',0.2);
% I = imnoise(G,'speckle',0.2);
I=double(I);


%%
v=[90  150];
c=2;
% v=[50 100 200];
% c=3;
% v=[40 90 150 255];
% c=4;
% v=[10 50 100 175 255];
% c=5;

m1=2;
% [~,~,tp2]=fcm(c,v,double(G),m1);
% figure,imshow(uint8(tp2))

%% fcm
[u,vi,I1]=fcm(c,v,I,m1);


%% KWFLICM_WGBF
tt=clock;
r=7;
im=GFF(I,3,5);
da=0.1; 
[outu10,~]=RDKWFLICM(im,I,c,m1,da,r,v);
img9=label1(outu10,c);
img9=img9(r+1:end-r,r+1:end-r);
figure;imshow(uint8(img9));title('KWFLICM_WGBF'); 
% disp(['KWFLICM_WGBF运行时间：',num2str(etime(clock,tt))]); 
disp(num2str(etime(clock,tt))); 
