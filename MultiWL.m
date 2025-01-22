function [fa,aq,ap,Tqp] = MultiWL(x,y,detaq,Q,P)
%% 计算输入的Leader
lenx = length(x);

[Lotmp,Hitmp] = wfilters('db3');       %获取对应小波的高低通滤波系数
Lo = coder.const(cast(Lotmp,'like',x));    %
Hi = coder.const(cast(Hitmp,'like',x));    %参数设置提高生成代码的效率和性能

Ntotal = lenx+length(Lotmp)-1;             %计算卷积总长度为dwtleaders函数提供输入

[allleadersx,~,ncount] = dwtleaders(x,coder.const(Lo),coder.const(Hi),Ntotal);  %通过dwtleaders函数获取输入x的leader
[allleadersy,~,~] = dwtleaders(y,coder.const(Lo),coder.const(Hi),Ntotal);       %通过dwtleaders函数获取输入y的leader

Nest = length(ncount);   %对尺度数量进行计数
%% 构造并计算 结构函数multivariate structure function Sf((q,p),j)
Sq = [];               %用于存放计算的结构函数
q = -Q:detaq:Q;        %
q = q(:);              %生成q,p序列，便于结构函数的构造
p = -P:detaq:P;        %
for j = 1:Nest
    numcoefsx = numel(allleadersx{j});  %当前尺度下的leader数量


    dxq = abs(allleadersx{j}).^q;            %通过广播机制生成xleader在尺度j下的所有q次幂的矩阵
    dyp = abs(allleadersy{j})'.^p;           %通过广播机制生成yleader在尺度j下的所有p次幂的矩阵

    zetaq = dxq*dyp;                    %通过矩阵乘法得到在对应（q,p）下的结构函数求和部分
    zetaq = zetaq/(numcoefsx);          %结构函数求和部分前的1/2^-j,替换为leader数目
    %zetaq = zetaq*2^(i);               %
 
    Sq = cat(3, Sq, zetaq);             %在第三维堆叠不同尺度下的结构函数
end

%% 线性拟合 缩放函数scaling function η(q,p)

Tqp = [];
xx = 1:Nest;
q_n = 1;
for x = -Q:detaq:Q
    Tp = [];
    p_n = 1;
    for y = -P:detaq:P
        yy = reshape(log2(Sq(q_n,p_n,:)),size(xx));  %根据论文中缩放函数的公式
        l = polyfit(xx,yy',1);                       %在尺度维上线性拟合得到对应缩放函数的数值
        Tp = [ Tp l(1)];                             % 进行线性拟合之前是否有什么办法可以预处理一下，去掉错误的点呢
        p_n = p_n+1;
    end
    q_n = q_n+1;
    Tqp = [Tqp;Tp];
end
%% Legendre transform

[fa,aq,ap] = mylegendreM(Tqp,detaq);    %legendre变换函数

end

