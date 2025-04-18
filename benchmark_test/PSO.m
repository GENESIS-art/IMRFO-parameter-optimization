function [BestX, BestF, HisBestF] = PSO(FunIndex, MaxIt, nPop)
    % 获取搜索空间范围和维度
    [Low, Up, Dim] = FunRange(FunIndex);

    % 初始化粒子位置和速度
    empty_particle.Position = [];
    empty_particle.Velocity = [];
    empty_particle.Cost = [];
    empty_particle.Best.Position = [];
    empty_particle.Best.Cost = inf;

    particle = repmat(empty_particle, nPop, 1);
    GlobalBest.Cost = inf;

    for i = 1:nPop
        % 初始化位置和速度
        particle(i).Position = rand(1, Dim) .* (Up - Low) + Low;
        particle(i).Velocity = zeros(1, Dim);

        % 计算适应度值
        particle(i).Cost = BenFunctions(particle(i).Position, FunIndex, Dim);

        % 更新个体最优解
        particle(i).Best.Position = particle(i).Position;
        particle(i).Best.Cost = particle(i).Cost;

        % 更新全局最优解
        if particle(i).Best.Cost < GlobalBest.Cost
            GlobalBest = particle(i).Best;
        end
    end

    % 记录历史最优适应度值
    HisBestF = zeros(MaxIt, 1);

    % PSO 参数
    w = 0.7; % 惯性权重
    c1 = 1.5; % 个体学习因子
    c2 = 1.5; % 社会学习因子

    % 迭代优化
    for it = 1:MaxIt
        for i = 1:nPop
            % 更新速度
            particle(i).Velocity = w * particle(i).Velocity + ...
                c1 * rand(1, Dim) .* (particle(i).Best.Position - particle(i).Position) + ...
                c2 * rand(1, Dim) .* (GlobalBest.Position - particle(i).Position);

            % 更新位置
            particle(i).Position = particle(i).Position + particle(i).Velocity;

            % 边界处理
            particle(i).Position = SpaceBound(particle(i).Position, Up, Low);

            % 计算适应度值
            particle(i).Cost = BenFunctions(particle(i).Position, FunIndex, Dim);

            % 更新个体最优解
            if particle(i).Cost < particle(i).Best.Cost
                particle(i).Best.Position = particle(i).Position;
                particle(i).Best.Cost = particle(i).Cost;

                % 更新全局最优解
                if particle(i).Best.Cost < GlobalBest.Cost
                    GlobalBest = particle(i).Best;
                end
            end
        end

        % 记录历史最优适应度值
        HisBestF(it) = GlobalBest.Cost;
    end

    % 返回结果
    BestX = GlobalBest.Position;
    BestF = GlobalBest.Cost;
end