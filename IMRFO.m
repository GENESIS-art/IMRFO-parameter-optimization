function [BestX,BestF,HisBestFit]=IMRFO(F_index,MaxIt,nPop)

         [Low,Up,Dim]=FunRange(F_index); 
iscontinue=0;
  if iscontinue==1
        load('BestResults.mat', 'PopPos', 'PopFit', 'BestX', 'BestF', 'BestX_history', 'BestF_history', 'HisBestFit','LastIteration');
            fprintf('Loaded saved data from previous run. Resuming from iteration %d.\n', LastIteration + 1);
            startIt = LastIteration + 1; 
  else
              
        if ~exist('BestX_history', 'var') || isempty(BestX_history)
            BestX_history = zeros(MaxIt, Dim);
        end
        if ~exist('BestF_history', 'var') || isempty(BestF_history)
            BestF_history = inf(MaxIt, 1);
        end
        if ~exist('HisBestFit', 'var') || isempty(HisBestFit)
            HisBestFit = zeros(MaxIt, 1);
        end   



        % 参数定义

        
        % 初始化
        Map = rand(nPop, Dim);  
            for d = 1:Dim  % 遍历所有维度
                a = 0.45+0.1*rand;  % 设定一个合适的a值0.45-0.55
                b = 0.2+0.02*rand; % 设定一个合适的b值
                for i = 2:nPop  % 遍历所有个体
                    % 使用圆形映射公式更新位置
                    Map(i, d) = mod(Map(i-1, d) + b - (a / (2 * pi) ) * sin(2 * pi * Map(i-1, d)), 1);
                end
            end

%         % 绘制 Map 的前两列
%         figure;
%         scatter(1:nPop, Map(:,1), 10,'filled');
%         xlabel('个体');
%         ylabel('Circle映射值');
% 
%         figure;
%         scatter(1:nPop, Map(:,2), 10,'filled');
%         xlabel('个体');
%         ylabel('Circle映射值');

        PopPos = Low + Map .* (Up - Low);

%         if Dim == 2
% 
%             figure;
%             scatter(PopPos(:, 1), PopPos(:, 2), 10, 'filled');
%             xlabel('维度1');
%             ylabel('维度2');
% 
%             grid on;
%         elseif Dim > 2
% 
%             figure;
%             scatter3(PopPos(:, 1), PopPos(:, 2), PopPos(:, 3), 50, 'filled');
%             xlabel('X');
%             ylabel('Y');
%             zlabel('Z');
%             title('3D Population Distribution Using Tent Mapping');
%             grid on;
%             axis equal;
%         end

        for i=1:nPop   
%             PopPos(i,:)=rand(1,Dim).*(Up-Low)+Low;
            PopFit(i)=BenFunctions(PopPos(i,:),F_index,Dim);   
        end
           BestF=inf;
           BestX=[];
    
        for i=1:nPop
            if PopFit(i)<=BestF
               BestF=PopFit(i);
               BestX=PopPos(i,:);
            end
        end
        startIt = 1; 
        LastIteration = 0;
        save('BestResults.mat', 'PopPos', 'PopFit', 'BestX', 'BestF', 'BestX_history', 'BestF_history', 'HisBestFit', 'LastIteration');

  end

for It=startIt:MaxIt  
     clc;
     Coef=It/MaxIt; 
     P = 0.5 * (1 / (1 + exp(-(12 * Coef - 6)))) + 0.25;
    
       if rand<P
          r1=rand;                         
          Beta=2*exp(r1*((MaxIt-It+1)/MaxIt))*(sin(2*pi*r1));    
          if  Coef>rand  
            
              newPopPos(1,:)=BestX+rand(1,Dim).*(BestX-PopPos(1,:))+Beta*(BestX-PopPos(1,:)); %Equation (4)
          else
        
              IndivRand=rand(1,Dim).*(Up-Low)+Low;                                
              newPopPos(1,:)=IndivRand+rand(1,Dim).*(IndivRand-PopPos(1,:))+Beta*(IndivRand-PopPos(1,:)); %Equation (7)         
          end              
       else 
 
            Alpha=2*rand(1,Dim).*(-log(rand(1,Dim))).^0.5;           
            newPopPos(1,:)=PopPos(1,:)+rand(1,Dim).*(BestX-PopPos(1,:))+Alpha.*(BestX-PopPos(1,:)); %Equation (1)
       end

    for i=2:nPop
        if rand<0.5
           r1=rand;                         
           Beta=2*exp(r1*((MaxIt-It+1)/MaxIt))*(sin(2*pi*r1));    
             if  Coef>rand                                                      
                 newPopPos(i,:)=BestX+rand(1,Dim).*(PopPos(i-1,:)-PopPos(i,:))+Beta*(BestX-PopPos(i,:)); %Equation (4)
             else
                 IndivRand=rand(1,Dim).*(Up-Low)+Low;                                
                 newPopPos(i,:)=IndivRand+rand(1,Dim).*(PopPos(i-1,:)-PopPos(i,:))+Beta*(IndivRand-PopPos(i,:));  %Equation (7)       
             end              
        else
            Alpha=2*rand(1,Dim).*(-log(rand(1,Dim))).^0.5;           
            newPopPos(i,:)=PopPos(i,:)+rand(1,Dim).*(PopPos(i-1,:)-PopPos(i,:))+Alpha.*(BestX-PopPos(i,:)); %Equation (1)
       end         
    end
         
           for i=1:nPop        
               newPopPos(i,:)=SpaceBound(newPopPos(i,:),Up,Low);
               newPopFit(i)=BenFunctions(newPopPos(i,:),F_index,Dim);    
              if newPopFit(i)<PopFit(i)
                 PopFit(i)=newPopFit(i);
                 PopPos(i,:)=newPopPos(i,:);
              end
           end



            S=2;
        for i=1:nPop
              levy_step = slight(PopPos(i, :));
              newPopPos(i,:)=PopPos(i,:)+(S).*(rand*BestX-rand*PopPos(i,:))+0.5*levy_step(1,:).*(BestX-PopPos(i,:)); 
        end
     for i=1:nPop        
         newPopPos(i,:)=SpaceBound(newPopPos(i,:),Up,Low);
         newPopFit(i)=BenFunctions(newPopPos(i,:),F_index,Dim);    
         if newPopFit(i)<PopFit(i)
            PopFit(i)=newPopFit(i);
            PopPos(i,:)=newPopPos(i,:);
         end
     end
     
     for i=1:nPop
        if PopFit(i)<BestF
           BestF=PopFit(i);
           BestX=PopPos(i,:);            
        end
     end

        BestX_history(It, :) = BestX;
        BestF_history(It) = BestF;
        HisBestFit(It) = BestF;
        LastIteration = It; 


%         save('BestResults.mat', 'PopPos', 'PopFit', 'BestX', 'BestF', 'BestX_history', 'BestF_history', 'HisBestFit', 'LastIteration', '-append');
 end
     HisBestFit(It)=BestF;
end



function step = slight(x)

    beta = 1.1; 
    sigma_u = (gamma(1 + beta) * sin(pi * beta / 2)) / (gamma((1 + beta) / 2) * beta * 2^(beta / 2));
    sigma_v = 1; 


    u = randn(size(x)) * sigma_u;
    v = randn(size(x)) * sigma_v;


    step = u ./ (abs(v).^(1 / beta));  
end
