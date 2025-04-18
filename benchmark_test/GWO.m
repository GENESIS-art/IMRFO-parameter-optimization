function [BestX, BestF, HisBestF] = GWO(FunIndex, MaxIt, nPop)
    % 获取搜索空间范围和维度
    [Low, Up, Dim] = FunRange(FunIndex);

    % 初始化狼群位置
    Alpha = zeros(1, Dim);
    Beta = zeros(1, Dim);
    Delta = zeros(1, Dim);
    AlphaScore = inf;
    BetaScore = inf;
    DeltaScore = inf;

    % 初始化种群位置
    PopPos = rand(nPop, Dim) .* (Up - Low) + Low;
    PopFit = zeros(nPop, 1);

    % 计算种群适应度
    for i = 1:nPop
        PopFit(i) = BenFunctions(PopPos(i, :), FunIndex, Dim);
    end

    % 记录历史最优适应度值
    HisBestF = zeros(MaxIt, 1);

    % 初始化全局最优解
    GlobalBestX = PopPos(1, :);  % 初始化全局最优位置
    GlobalBestF = PopFit(1);      % 初始化全局最优适应度

    % 迭代优化
    for t = 1:MaxIt
        % 更新Alpha、Beta、Delta
        [~, SortIndex] = sort(PopFit);
        Alpha = PopPos(SortIndex(1), :);
        AlphaScore = PopFit(SortIndex(1));
        Beta = PopPos(SortIndex(2), :);
        BetaScore = PopFit(SortIndex(2));
        Delta = PopPos(SortIndex(3), :);
        DeltaScore = PopFit(SortIndex(3));

        % 更新种群位置
        a = 2 - t * (2 / MaxIt); % 线性递减的a值
        for i = 1:nPop
            A1 = 2 * a * rand(1, Dim) - a;
            C1 = 2 * rand(1, Dim);
            DAlpha = abs(C1 .* Alpha - PopPos(i, :));
            X1 = Alpha - A1 .* DAlpha;

            A2 = 2 * a * rand(1, Dim) - a;
            C2 = 2 * rand(1, Dim);
            DBeta = abs(C2 .* Beta - PopPos(i, :));
            X2 = Beta - A2 .* DBeta;

            A3 = 2 * a * rand(1, Dim) - a;
            C3 = 2 * rand(1, Dim);
            DDelta = abs(C3 .* Delta - PopPos(i, :));
            X3 = Delta - A3 .* DDelta;

            PopPos(i, :) = (X1 + X2 + X3) / 3;

            % 边界处理
            PopPos(i, :) = SpaceBound(PopPos(i, :), Up, Low);

            % 计算适应度值
            PopFit(i) = BenFunctions(PopPos(i, :), FunIndex, Dim);
        end

        % 获取当前最优适应度和位置
        [BestF, idx] = min(PopFit);
        BestX = PopPos(idx, :);

        % 更新全局最优解
        if BestF < GlobalBestF
            GlobalBestF = BestF;
            GlobalBestX = BestX;
        end

        % 记录历史最优适应度值
        if t == 1
            % 第一次迭代，直接赋值
            HisBestF(t) = BestF;
        else
            % 不是第一次迭代，比较当前最优值与历史最优值
            if BestF < min(HisBestF(1:t-1))
                HisBestF(t) = BestF;
            else
                HisBestF(t) = min(HisBestF(1:t-1));
            end
        end
    end

    % 最终全局最优解
    BestX = GlobalBestX;
    BestF = GlobalBestF;
end
