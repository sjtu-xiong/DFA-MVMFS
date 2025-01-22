function [Lh,hq,hp,tq] = DFAMFS(signal1,signal2,scale,m,q)
% signal:       input signal
% scale:        vector of scales
% q:            q-order that weights the local variations 
% m:            polynomial order for the detrending
%
% Hq:           q-order Hurst exponent  Hurst 指数
% tq:           q-order mass exponent   质量指数
% hq:           q-order singularity exponent   奇异指数
% Dq:           q-order dimension   维数
% Fq:           q-order scaling function    尺度函数

% EXAMPLE------------------------------------------------------------------
% scmin=10;
% scmax=410;
% scres=10;
% exponents=linspace(log2(scmin),log2(scmax),scres);
% scale=round(2.^exponents);
% q=linspace(-5,5,41);
% m=1;

%% 对输入处理
warning off
X=cumsum(signal1-mean(signal1)); %构造新的x(t)序列

Y=cumsum(signal2-mean(signal2)); %构造新的y(t)序列
p=q;
if min(size(X))~=1||min(size(scale))~=1||min(size(q))~=1
    error('Input arguments signal, scale and q must be a vector');
end
if size(X,2)==1
   X=transpose(X);
end

if size(Y,2)==1
   Y=transpose(Y);
end

%检查是否为行向量，并将其设置为行向量
if min(scale)<m+1
   error('The minimum scale must be larger than trend order m+1')
end
%% 序列构造 求解缩放函数scaling function T(q,p)
for ns=1:length(scale)
    segments(ns)=floor(length(X)/scale(ns));
    for v=1:segments(ns)
        Index=((((v-1)*scale(ns))+1):(v*scale(ns)));
        CX=polyfit(Index,X(Index),m);      %m阶多项式拟合
        CY=polyfit(Index,Y(Index),m);      %
        fitX=polyval(CX,Index);            %
        fitY=polyval(CY,Index);            %
        RMS_scaleX{ns}(v)=sqrt(mean((X(Index)-fitX).^2)); %Fx(v,s)
        RMS_scaleY{ns}(v)=sqrt(mean((Y(Index)-fitY).^2)); %Fy(v,s)
    end
    for nq=1:length(q)
        for np=1:length(p)
            qRMS{nq,np,ns}=(RMS_scaleX{ns}.^q(nq)).*(RMS_scaleY{ns}.^p(np));   % Fq,p(s)
            Fq(nq,np,ns) = sum(qRMS{nq,np,ns});

            %Fq(nq,ns)=mean(qRMS{nq,ns}).^(1/q(nq));
            % Fq(nq,ns)=mean(qRMS{nq,ns});
        end
    end
end
for nq=1:length(q)
    for np=1:length(p)
        C = polyfit(log2(scale),log2(Fq(nq,np,:)),1); %尺度上拟合
        tq(nq,np) = C(1);
    end
end
%% Legendre transform
[Lh,hq,hp] = mylegendreM(tq,q(2)-q(1));