function [outu,label]=KWFLICM_WGBF_RGB(G,I0,c,m1,da,r,v)
im=double(G);
%%参数
ee=0;
kk=0;
dr=1;
T=4;
%%
sigma=1000;
I0=GFF(I0,3,0.05);
I0=uint8(I0);
I0=padarray(I0,[r r 0],'replicate');

if ~isa(I0,'double')
    I0 = double(I0)/255;   % 转换为double类型
end
[m,n,q]=size(I0);
im=uint8(im);
P=padarray(im,[r r 0],'replicate');
if ~isa(im,'double')
 P = double(P)/255;
end
%%
f=zeros(m,n,q);
gk=GaussianKernel(r,2);
v=v./255;
% [u,vi,~]=fcm(c,v,I0,m1);
[vi,u,~,~]=fcm_RGB(c,I0,m1,v);
% fcm_RGB(c,I,m1,v);
%%gr wsc
gr=zeros(2*r+1);
% wsc=zeros(2*r+1);
wsc=zeros(3);
for ii=1:2*r+1
    for jj=1:2*r+1
        gr(ii,jj)=exp(-((ii-r-1)^2+(jj-r-1)^2)/(2*r^2)); % 欧氏距离
        
    end
end
tp=zeros(m,n);
var=zeros(m,n);
%%初始化fij
for i=(r+1):m-r
    for j=(r+1):n-r
        x=-r:1:r;
        y=-r:1:r;
        %range kernel
        tin=I0(i+x,j+y,:);       
        win1=P(i+x,j+x,:);% 获取窗口
        win1=(win1(r+1,r+1,:)-win1); % 相对中心点处的强度差异，中心点为参考灰度值
        gs1=exp(-(win1.^2)/(2*(da^2))); % 高斯函数
        for k=1:q
        weights(:,:,k)=(gr.*gs1(:,:,k))/sum(sum(gr.*gs1(:,:,k))); % 高斯权值
        end
        f(i,j,:) =sum(sum(weights.*tin));
    end
end 
while ee<0.0001&&kk<T
    v=vi;
    
    gk1=GK(f,v,sigma,c);%%高斯核
    %%更新fij
    for i=(r+1):m-r
        for j=(r+1):n-r
            x=-r:1:r;
            y=-r:1:r;
            w=sum(abs(P(i+x,j+y,:)-P(i,j,:)))./3;
            gs=exp(-(w.^2)/(2*(da^2))); % 高斯函数
            %%fi
            t1=0;t2=0;
            for k=1:c
                t1=t1+u(i,j,k)^m1.*v(:,k).*gk1(i,j,k);
                t2=t2+u(i,j,k)^m1*gk1(i,j,k);
            end
            ttpp=sum(gr.*gs.*I0(i+x,j+y,:),1);
            f(i,j,:)=(t1+reshape(dr.*sum(ttpp,2),[3,1]))./  ...
                   (t2+dr.*sum(sum(gr.*reshape(gs,[3 3]))));
        end
    end

    gk=GK(f,v,sigma,c);%%高斯核

    for i=(r+1):m-r
        for j=(r+1):n-r
            x=-1:1;
            y=-1:1:1;
            %%var
            tp1=sum(sum(sum(f(i+x,j+y,:))))./3;
            t1=tp1/(9);%%邻域均值
            var(i,j)=sum(sum((sum(f(i+x,j+y,:),3)./3-t1).^2))/9;%%方差
            tp(i,j)=var(i,j)/t1^2;%%cj
        end
    end
    for i=(r+1):m-r
        for j=(r+1):n-r
            x=-1:1:1;
            y=-1:1:1;
            %%ksi
            t2=sum(sum(tp(i+x,j+y)))/((2*r+1)^2);%%c-
            ksi=exp(-(tp(i+x,j+y)-t2));%%r*r
            eta=ksi./sum(sum(ksi));%%
            wgc=zeros(3);
            for ii=-1:1
                for jj=-1:1
                    wsc(ii+2,jj+2)=1/(1+ii^2+jj^2);
                    if tp(ii+i,jj+j)<t2
                        wgc(ii+2,jj+2)=2+eta(ii+2,jj+2);
                    else
                        wgc(ii+2,jj+2)=2-eta(ii+2,jj+2);
                    end
                end
            end            
            wij=wsc.*wgc;
            %%uij
            t3=wij.*((1-u(i+x,j+y,:)).^m1).*(1-gk(i+x,j+y,:));
            t4=(1-gk(i,j,:)+sum(sum(t3))).^(-1/(m1-1));
            u(i,j,:)=t4/(sum(t4)+0.0001);
%              u(i,j,k)=1/((f(i,j)-v(k))^2+0.0001+G(i,j,k))^(1/(m1-1))*1/t;
        end
    end
    %%vk
    for k=1:c
        tp2=0.0;
        tp3=0.0;
        for i=(r+1):m-r
            for j=(r+1):n-r
                tp2=tp2+u(i,j,k).^m1*gk(i,j,k)*f(i,j,:);
                tp3=tp3+u(i,j,k).^m1*gk(i,j,k);
            end
        end 
        vi(:,k)=tp2./(tp3+0.0001);
    end
    %%终止条件
    temp=0.0;
    for k=1:c
        temp=temp+(v(:,k)-vi(:,k)).^2;
    end
    if   temp < 0.0001
        ee=0.0001;
    end
    kk=kk+1;
%  UUUUU=u(79:79+2,191:191+2,:)
% FFFFF=f(79:79+2,191:191+2)
% VVVV=vi*255


end



outu=u;



label=zeros(m,n);

for i=1:m
    for j=1:n        
        if u(i,j,1)>u(i,j,2) && u(i,j,1)>u(i,j,3)
            label(i,j)=1;
        elseif u(i,j,2)>u(i,j,1) && u(i,j,2)>u(i,j,3)
            label(i,j)=2;
        elseif u(i,j,3)>u(i,j,1) && u(i,j,3)>u(i,j,2)
            label(i,j)=3;
        end        
    end
end


end
 

function I=GK(x1,x2,sigma,c)
[m,n,~]=size(x1);
% I=zeros(3,3);
% d=zeros(3,3);
for i=1:m
    for j=1:n
        for k=1:c
        d(i,j,k)=sum((reshape(x1(i,j,:),[1,3])-x2(:,k)').^2,2)/3+0.0001;
        I(i,j,k)=exp(-d(i,j,k)/(2*sigma^2));        
        end
    end
end

end