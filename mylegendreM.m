function [Lh,hq,hp] = mylegendreM(Tq,detaq)
%% paprameter descri.
% legendre transform
% Tq is the mass function, q is the moment and aq is the singularity
% exponent coming from the derivation of q. and fa is the MFS: fa =
% q*a(q)-Tq;
% Ĭ��q�ļ����1��Ȼ��qĬ���Ǵ�-Q:1��Q�������ڲ�ֵ�ʱ��ֱ��ѡ�����2��
%% null variable.
hq = [];
hp = [];
%% ����������ָ��  H(q) =  dTq/dq
[Q,P] = size(Tq);
for j = 1:1:P
    hq1=[];
    for i = 2:1:Q-1
        hq1 = [hq1;(Tq(i+1,j)-Tq(i-1,j))/2];
    end
    hq1 = [(Tq(2,j)-Tq(1,j));hq1;(Tq(Q,j)-Tq(Q-1,j))];
    hq = [hq,hq1];
end

for j = 1:1:Q
    hp1=[];
    for i = 2:1:P-1
        hp1 = [hp1,(Tq(j,i+1)-Tq(j,i-1))/2];
    end
    hp1 = [(Tq(j,2)-Tq(j,1)),hp1,(Tq(j,Q)-Tq(j,Q-1))];
    hp = [hp;hp1];
end
hq = hq./(detaq);
hp = hp./(detaq);
%% ������ط�����Lh=1+hq*q+hp*p-Tq
q = (-(Q-1)/2):1:(Q-1)/2;    
q = detaq.*q;
p =(-(P-1)/2):1:(P-1)/2; 
p = detaq.*p;
Lh = hq.*q'+hp.*p-Tq;      %DFA-MFS
% Lh = 1+hq.*q'+hp.*p-Tq;      %WL-MFS
Lh = real(Lh); 

