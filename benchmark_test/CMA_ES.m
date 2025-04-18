function [BestX, BestF, HisBestFit] = CMAES(F_index, MaxIt, nPop)
    [Low, Up, Dim] = FunRange(F_index);
    sigma = 0.3 * (Up - Low);  % 初始步长，设为搜索空间范围的30%
    mu = nPop / 2;             % 父代数量
    lambda = nPop;             % 子代数量
    weights = log(mu + 0.5) - log(1:mu)';
    weights = weights / sum(weights);  % 权重归一化
    mu_eff = 1 / sum(weights.^2);      % 有效父代数

    % 初始化策略参数
    cs = (mu_eff + 2) / (Dim + mu_eff + 5);
    ds = 1 + cs + 2 * max(sqrt((mu_eff - 1) / (Dim + 1)) - 1, 0);
    cc = (4 + mu_eff / Dim) / (Dim + 4 + 2 * mu_eff / Dim);
    c1 = 2 / ((Dim + 1.3)^2 + mu_eff);
    cmu = min(1 - c1, 2 * (mu_eff - 2 + 1 / mu_eff) / ((Dim + 2)^2 + mu_eff));
    chiN = sqrt(Dim) * (1 - 1 / (4 * Dim) + 1 / (21 * Dim^2));

    % 初始化动态参数
    m = Low + rand(1, Dim) .* (Up - Low);  % 初始均值
    C = eye(Dim);                          % 协方差矩阵
    epsilon = 1e-8;
    C = C + epsilon * eye(Dim);
    C = (C + C') / 2;
    pc = zeros(1, Dim);                    % 路径向量
    ps = zeros(1, Dim);                    % 步长路径
    HisBestFit = zeros(MaxIt, 1);          % 最佳适应度历史
    BestF = inf;
    BestX = m;

    for It = 1:MaxIt
        % 采样子代
        arz = randn(lambda, Dim);
        [U, S, ~] = svd(C);
        ary = arz * U * sqrt(S);
        arx = m + sigma .* ary;
        % 边界处理
        arx = min(max(arx, Low), Up);
        % 评估适应度
        fitness = zeros(lambda, 1);
        for k = 1:lambda
            fitness(k) = BenFunctions(arx(k, :), F_index, Dim);
        end
        % 选择前 mu 个最优个体
        [~, idx] = sort(fitness);
        x_selected = arx(idx(1:mu), :);
        y_selected = ary(idx(1:mu), :);
        % 更新均值
        m_old = m;
        m = weights' * x_selected;
        % 更新路径
        y_w = weights' * y_selected;
        ps = (1 - cs) * ps + sqrt(cs * (2 - cs) * mu_eff) * y_w / sigma;
        hsig = norm(ps) / sqrt(1 - (1 - cs)^(2 * It)) / chiN < 1.4 + 2 / (Dim + 1);
        pc = (1 - cc) * pc + hsig * sqrt(cc * (2 - cc) * mu_eff) * y_w;
        % 更新协方差矩阵
        C = (1 - c1 - cmu) * C + ...
            c1 * (pc' * pc + (1 - hsig) * cc * (2 - cc) * C) + ...
            cmu * (y_selected' * diag(weights) * y_selected);
        % 更新步长
        sigma = sigma * exp((cs / ds) * (norm(ps) / chiN - 1));
        % 更新历史最佳
        if fitness(idx(1)) < BestF
            BestF = fitness(idx(1));
            BestX = arx(idx(1), :);
        end
        HisBestFit(It) = BestF;
    end
end  

