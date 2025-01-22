% 参数设定
p1 = 0.1;  % 第一个信号的概率 p1
p2 = 0.4;  % 第二个信号的概率 p2
depth = 12; % 迭代的层数

% 生成两个随机选择的二项测度信号
signal1 = binomial_measure(p1, depth);
signal2 = binomial_measure(p2, depth);

% 绘信号
x1 = linspace(0, 1, length(signal1));  % 生成x轴坐标 for signal1
x2 = linspace(0, 1, length(signal2));  % 生成x轴坐标 for signal2

% 使用subplot绘制两个信号，分别显示在上下部分
figure;

% 第一个信号
subplot(2,1,1);
plot(x1, signal1, 'LineWidth', 1.5);
title('Randomized Binomial Measure Signal 1');
xlabel('x');
ylabel('Measure');

% 第二个信号
subplot(2,1,2);
plot(x2, signal2, 'LineWidth', 1.5);
title('Randomized Binomial Measure Signal 2');
xlabel('x');
ylabel('Measure');
% 
% 函数定义
% 
function measure = binomial_measure(p, depth)
% 生成随机权重选择的二项测度信号
% 参数:
%   p - 给定的概率
%   depth - 迭代层数
% 返回:
%   生成的测度信号数组
%
% 初始测度
measure = 1;

% 递归迭代生成二项测度
for i = 1:depth
    %     随机决定左分支为 p 还是 1-p
    if rand < 0.5
        left_weight = p;
        right_weight = 1 - p;
    else
        left_weight = 1 - p;
        right_weight = p;
    end

    % 按照随机权重生成新一层的测度
    measure = [left_weight * measure, right_weight * measure];
end
end


