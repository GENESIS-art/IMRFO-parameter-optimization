function [BestX, BestF, HisBestFit] = ALO(F_index, MaxIt, nPop)
    % 输入参数:
    % F_index: 测试函数编号
    % MaxIt: 最大迭代次数
    % nPop: 种群大小
    %
    % 输出参数:
    % BestX: 最优解的位置
    % BestF: 最优解的目标函数值
    % HisBestFit: 每次迭代的最优适应度历史记录

    % 获取问题的边界和维度
    [Low, Up, Dim] = FunRange(F_index);
    
    % 初始化种群
    antlion = zeros(nPop, Dim);   % 蚁狮位置
    ant = zeros(nPop, Dim);       % 蚂蚁位置
    antlion_fit = zeros(nPop, 1); % 蚁狮适应度
    ant_fit = zeros(nPop, 1);     % 蚂蚁适应度
    
    % 初始化蚁狮和蚂蚁的位置
    for i = 1:nPop
        antlion(i,:) = rand(1,Dim).*(Up-Low) + Low;
        ant(i,:) = rand(1,Dim).*(Up-Low) + Low;
        
        % 计算初始适应度
        antlion_fit(i) = BenFunctions(antlion(i,:), F_index, Dim);
        ant_fit(i) = BenFunctions(ant(i,:), F_index, Dim);
    end
    
    % 排序并找到最优蚁狮
    [antlion_fit, sort_index] = sort(antlion_fit);
    antlion = antlion(sort_index, :);
    
    BestF = antlion_fit(1);
    BestX = antlion(1,:);
    
    HisBestFit = zeros(MaxIt, 1);
    
    % 开始迭代
    for It = 1:MaxIt
        % 精英保留策略
        elite_antlion = antlion(1,:);
        elite_fit = antlion_fit(1);
        
        % 蚂蚁围绕蚁狮随机游走
        for i = 1:nPop
            % 选择蚁狮(轮盘赌选择)
            Rolette_index = RouletteWheelSelection(1./antlion_fit);
            if Rolette_index == -1
                Rolette_index = 1;
            end
            
            % 计算游走的边界
            I = 1; % 缩放因子
            if It > 0.1*MaxIt
                I = 1 + 100*(It/MaxIt);
            end
            if It > 0.5*MaxIt
                I = 1 + 1000*(It/MaxIt);
            end
            if It > 0.75*MaxIt
                I = 1 + 10000*(It/MaxIt);
            end
            if It > 0.9*MaxIt
                I = 1 + 100000*(It/MaxIt);
            end
            if It > 0.95*MaxIt
                I = 1 + 1000000*(It/MaxIt);
            end
            
            % 更新边界
            lb = antlion(Rolette_index,:) ./ I;
            ub = antlion(Rolette_index,:) .* I;
            
            % 确保边界不超出问题边界
            lb = max(lb, Low);
            ub = min(ub, Up);
            
            % 随机游走
            if rand < 0.5
                ant(i,:) = (antlion(Rolette_index,:) + lb) .* rand(1,Dim) + (antlion(Rolette_index,:) - lb) .* rand(1,Dim);
            else
                ant(i,:) = (antlion(Rolette_index,:) + ub) .* rand(1,Dim) + (antlion(Rolette_index,:) - ub) .* rand(1,Dim);
            end
            
            % 边界检查
            ant(i,:) = SpaceBound(ant(i,:), Up, Low);
            
            % 计算适应度
            ant_fit(i) = BenFunctions(ant(i,:), F_index, Dim);
        end
        
        % 蚂蚁与蚁狮竞争
        double_pop = [antlion; ant];
        double_fit = [antlion_fit; ant_fit];
        
        [double_fit, sort_index] = sort(double_fit);
        double_pop = double_pop(sort_index, :);
        
        % 选择前nPop个作为新的蚁狮
        antlion = double_pop(1:nPop, :);
        antlion_fit = double_fit(1:nPop);
        
        % 更新最优解
        if antlion_fit(1) < BestF
            BestF = antlion_fit(1);
            BestX = antlion(1,:);
        end
        
        % 记录历史最优适应度
        HisBestFit(It) = BestF;
        
        % 显示迭代信息
%         if mod(It, 50) == 0
%             disp(['Iteration ' num2str(It) ': Best Fitness = ' num2str(BestF)]);
%         end
    end
end

% 轮盘赌选择函数
function index = RouletteWheelSelection(weights)
    accumulation = cumsum(weights);
    p = rand() * accumulation(end);
    index = -1;
    for i = 1:length(accumulation)
        if accumulation(i) > p
            index = i;
            break;
        end
    end
end