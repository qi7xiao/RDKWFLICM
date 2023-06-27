function y=GFF(im,sigd,sigr)
% bilateral filter˫���˲���
% �������룺
%           im    �����ͼ��
%           sigd  �ռ��ں˵�ʱ�����
%           sigr  �ں˲���ǿ�ȱ仯��Χ
% ���������
%          out  �˲�ͼ�� = output imagespatial kernel
% sigr=(n*100)^2/(.003*(sigd^2));  % ����ӦRֵ��nΪ��˹����ǿ��,n=0.001
[~,~,channel]=size(im);
if channel==1
    y=calGFF(im,sigd,sigr);
elseif channel==3
   y(:,:,1)=calGFF(im(:,:,1),sigd,sigr);
   y(:,:,2)=calGFF(im(:,:,2),sigd,sigr);
   y(:,:,3)=calGFF(im(:,:,3),sigd,sigr);
end

function out=calGFF(im,sigd,sigr)
% ��˹Ƶ���˲��� 
        
% Padding ��չͼ��ı߽磬��ֹ�������ڱ߽�ֵ���
proci=padarray(im,[sigd sigd],'replicate');
[row,clm]=size(proci);    % Size of image
    proci = double(proci)/255;   % ת��Ϊdouble����
K=sigd;
L=-K:K;
c=K+1;   % ����Ԫ��λ��
for r=(1+K):(row-K)          % ��    
    for s=(1+K):(clm-K)      % ��                     
            win=proci((r+L),(s+L));% ��ȡ����   
            I=win; % �ҶȾ���
            win=win(c,c)-win; % ������ĵ㴦��ǿ�Ȳ��죬���ĵ�Ϊ�ο��Ҷ�ֵ
            win=sqrt(win.^2); % ��֤win�е�ÿһ��Ԫ��Ϊ��
            Gwi=exp(-(win.^2)/(2*(sigr^2))); % ��˹����      
            weights=(Gwi)/sum(sum(Gwi)); % ��˹Ȩֵ
            proci(r,s) =sum(sum(weights.*I));               % �õ���ǰ˫���˲�ֵ  
    end
end

% �Ƴ��߽���չֵ
out=im2uint8(proci(K+1:end-K,K+1:end-K));   % ����ת��


