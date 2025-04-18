%--------------------------------------------------------------------------
% WOA code v1.0.
% Developed in MATLAB R2011b
% The code is based on the following papers.
% S. Mirjalili, A. Lewis, The Whale Optimization Algorithm,
% Advances in Engineering Software, Volume 95, 2016, Pages 51-67,
% https://doi.org/10.1016/j.advengsoft.2016.01.008.
% -------------------------------------------------------------------------
% FunIndex: Index of function.
% MaxIt: The maximum number of iterations.
% PopSize: The size of population.
% PopPos: The position of population.
% PopFit: The fitness of population.
% Dim: The dimensionality of problem.
% a: The parameter that controls the spiral shape.
% BestF: The best fitness corresponding to BestX. 
% HisBestFit: History best fitness over iterations. 
% Low: The low bound of search space.
% Up: The up bound of search space.

function [BestX,BestF,HisBestFit]=WOA(F_index,MaxIt,nPop)

    [Low,Up,Dim]=FunRange(F_index); 

    % Initialize the positions of search agents
    for i=1:nPop   
        PopPos(i,:)=rand(1,Dim).*(Up-Low)+Low;
        PopFit(i)=BenFunctions(PopPos(i,:),F_index,Dim);   
    end

    % Initialize the best solution
    BestF=inf;
    BestX=[];

    for i=1:nPop
        if PopFit(i)<=BestF
           BestF=PopFit(i);
           BestX=PopPos(i,:);
        end
    end

    HisBestFit=zeros(MaxIt,1);

    % Main loop
    for It=1:MaxIt  
        % Update a, linearly decreases from 2 to 0
        a = 2 - It * (2 / MaxIt);
        
        for i=1:nPop
            % Update A, C, l, and p
            A = 2 * a * rand() - a;  % Equation (2.3)
            C = 2 * rand();          % Equation (2.4)
            l = (a-1)*rand()+1;       % Equation (2.5)
            p = rand();              % Probability of choosing encircling or spiral
            
            % Update the position of the current search agent
            if p < 0.5
                if abs(A) < 1
                    % Encircling prey (Equation 2.1)
                    D = abs(C .* BestX - PopPos(i,:));
                    newPopPos(i,:) = BestX - A .* D;
                else
                    % Search for prey (Equation 2.8)
                    randInd = randi([1, nPop]);
                    D = abs(C .* PopPos(randInd,:) - PopPos(i,:));
                    newPopPos(i,:) = PopPos(randInd,:) - A .* D;
                end
            else
                % Spiral updating position (Equation 2.6)
                D = abs(BestX - PopPos(i,:));
                newPopPos(i,:) = D .* exp(l) .* cos(2*pi*l) + BestX;
            end
            
            % Boundary check
            newPopPos(i,:)=SpaceBound(newPopPos(i,:),Up,Low);
            
            % Evaluate the new position
            newPopFit(i)=BenFunctions(newPopPos(i,:),F_index,Dim);    
            
            % Update the population if the new position is better
            if newPopFit(i)<PopFit(i)
                PopFit(i)=newPopFit(i);
                PopPos(i,:)=newPopPos(i,:);
            end
        end
        
        % Update the best solution
        for i=1:nPop
            if PopFit(i)<BestF
                BestF=PopFit(i);
                BestX=PopPos(i,:);            
            end
        end
        
        % Store the best fitness value
        HisBestFit(It)=BestF;
    end
end