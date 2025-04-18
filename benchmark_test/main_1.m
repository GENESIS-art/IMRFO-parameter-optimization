clc;
clear;

% 设置参数
MaxIteration = 2000 ; % 每轮最大迭代次数
PopSize = 100; % 种群大小
NumRounds = 10; % 循环10轮
FunIndex = 8; % 使用的目标函数
%1-F2;2-F5;3-F7;4-F8;
                %5-F9;6-F11;7-F12 ，10000轮迭代改进算法表现很好；2 5  8 12 13
% 初始化存储最佳值、平均值、最优值和标准差
BestValues = zeros(NumRounds, 7); % 存储每轮各个算法的最佳值
AllFitness = zeros(NumRounds, MaxIteration, 7); % 存储所有轮次的收敛曲线

% 运行算法
for round = 1:NumRounds
    % 运行不同算法
    [BestX1, BestF1, HisBestF1] = MRFO(FunIndex, MaxIteration, PopSize); % MRFO 算法
    [BestX2, BestF2, HisBestF2] = IIMRFO(FunIndex, MaxIteration, PopSize); % IMRFO 算法
    [BestX3, BestF3, HisBestF3] = PSO(FunIndex, MaxIteration, PopSize);  % PSO 算法
    [BestX4, BestF4, HisBestF4] = GWO(FunIndex, MaxIteration, PopSize);  % GWO 算法
    [BestX5, BestF5, HisBestF5] = WOA(FunIndex, MaxIteration, PopSize);   % WOA 算法
    [BestX6, BestF6, HisBestF6] = ALO(FunIndex, MaxIteration, PopSize);   % ALO 算法
    [BestX7, BestF7, HisBestF7] = CMA_ES(FunIndex, MaxIteration, PopSize); % CMA-ES算法
    

    % 存储收敛曲线和最佳值
    AllFitness(round, :, 1) = HisBestF1;
    AllFitness(round, :, 2) = HisBestF2;
    AllFitness(round, :, 3) = HisBestF3;
    AllFitness(round, :, 4) = HisBestF4;
    AllFitness(round, :, 5) = HisBestF5;
    AllFitness(round, :, 6) = HisBestF6;
    AllFitness(round, :, 7) = HisBestF7;

    
    BestValues(round, :) = [BestF1, BestF2, BestF3, BestF4, BestF5, BestF6, BestF7];

    % 输出当前轮次进度
    disp(['Round ', num2str(round), ' completed']);
end

% 计算IMRFO的最佳一轮（最佳适应度的那一轮）
[~, bestRoundMRFO] = min(BestValues(:, 1));
[~, bestRoundIMRFO] = min(BestValues(:, 2));
[~, bestRoundPSO] = min(BestValues(:, 3));
[~, bestRoundGWO] = min(BestValues(:, 4));
[~, bestRoundWOA] = min(BestValues(:, 5));
[~, bestRoundALO] = min(BestValues(:, 6));
[~, bestRoundCMAES] = min(BestValues(:, 7));
% 提取IMRFO最佳一轮的收敛曲线
BestF1_smooth = smoothdata(AllFitness(bestRoundMRFO, :, 1), 'movmean', 10);
BestF2_smooth = smoothdata(AllFitness(bestRoundIMRFO, :, 2), 'movmean', 10);
BestF3_smooth = smoothdata(AllFitness(bestRoundPSO, :, 3), 'movmean', 10);
BestF4_smooth = smoothdata(AllFitness(bestRoundGWO, :, 4), 'movmean', 10);
BestF5_smooth = smoothdata(AllFitness(bestRoundWOA, :, 5), 'movmean', 10);
BestF6_smooth = smoothdata(AllFitness(bestRoundALO, :, 6), 'movmean', 10);
BestF7_smooth = smoothdata(AllFitness(bestRoundCMAES, :, 7), 'movmean', 10);

% 绘制所有算法在IMRFO最佳一轮的收敛曲线
figure;
hold on;
semilogy(BestF1_smooth, 'g', 'LineWidth', 1.5); % MRFO 算法
semilogy(BestF2_smooth, 'r', 'LineWidth', 1.5); % IMRFO 算法
semilogy(BestF3_smooth, 'b', 'LineWidth', 1.5); % PSO 算法
semilogy(BestF4_smooth, 'm', 'LineWidth', 1.5); % GWO 算法
semilogy(BestF5_smooth, 'c', 'LineWidth', 1.5); % WOA 算法
semilogy(BestF6_smooth, 'y', 'LineWidth', 1.5); % ALO 算法
semilogy(BestF7_smooth, 'k', 'LineWidth', 1.5); % CMA-ES

hold off;

legend('MRFO', 'IMRFO', 'PSO', 'GWO', 'WOA', 'ALO', 'CMA-ES', 'FontSize', 10);
xlabel('迭代次数', 'FontSize', 14);  % 设置x轴标签的字体大小
ylabel('适应度值', 'FontSize', 14);  % 设置y轴标签的字体大小
% title(['F', num2str(FunIndex), ' - Smoothed Fitness Convergence'], 'FontSize', 16);  % 设置标题的字体大小
grid on;
set(gca, 'yscale', 'log');  % 设置y轴为对数尺度

% 计算平均值、最优值和标准差
MeanValues = mean(BestValues); % 平均值
MinValues = min(BestValues); % 最优值
StdValues = std(BestValues); % 标准差

% 输出平均值、最优值和标准差
disp('Mean values of best fitness for each algorithm:');
disp(MeanValues);
disp('Best values of fitness for each algorithm:');
disp(MinValues);
disp('Standard deviation of fitness for each algorithm:');
disp(StdValues);

% 绘制平均值柱状图
figure;
bar(MeanValues);
xticklabels({'MRFO', 'IMRFO', 'PSO', 'GWO', 'WOA', 'ALO', 'CMA-ES'});
ylabel('Average Fitness');
title('Average Fitness Comparison');
grid on;

% Wilcoxon秩和检验
alpha = 0.05;
IMRFO_vs_others = strings(1, 6);  % 用于记录与IMRFO对比结果的符号
p_values = zeros(1, 6);           % 保存每个比较的p值
z_values = zeros(1, 6);           % 保存signed rank值
% 提取每轮IMRFO的最优值
IMRFO_results = BestValues(:, 2);

% 与每个其他算法比较（1-MRFO，3-PSO，4-GWO，5-WOA，6-ALO）
alg_labels = {'MRFO', 'PSO', 'GWO', 'WOA', 'ALO', 'CMA-ES'};
for idx = 1:6
    alg_col = idx + (idx >= 2);  % 因为第2列是IMRFO，要跳过它
    other_results = BestValues(:, alg_col);
    [p, h, stats] = signrank(IMRFO_results, other_results);  % Wilcoxon符号秩检验
    p_values(idx) = p;  % 保存p值
    z_values(idx) = stats.signedrank;
    if p < alpha
        if median(IMRFO_results) < median(other_results)
            IMRFO_vs_others(idx) = "+";  % IMRFO显著更优
        else
            IMRFO_vs_others(idx) = "−";  % IMRFO显著更差
        end
    else
        IMRFO_vs_others(idx) = "=";      % 无显著差异
    end

    % 正确拼接字符串，避免列式输出
    disp(['IMRFO vs ', alg_labels{idx}, ': ', ...
          char(IMRFO_vs_others(idx)), ...
          ' (p = ', num2str(p, '%.4f'), ...
          ', signed rank = ', num2str(stats.signedrank), ')']);
end

% 最终结果表格输出
disp('Wilcoxon检验比较结果（IMRFO vs Others）:');
disp(table(alg_labels', IMRFO_vs_others', 'VariableNames', {'Algorithm', 'Result'}));
