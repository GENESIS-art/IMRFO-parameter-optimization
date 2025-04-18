
function [BestX,BestF,HisBestFit]=MRFO(F_index,MaxIt,nPop)

         [Low,Up,Dim]=FunRange(F_index); 

    for i=1:nPop   
        PopPos(i,:)=rand(1,Dim).*(Up-Low)+Low;
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

       HisBestFit=zeros(MaxIt,1);


for It=1:MaxIt  
     Coef=It/MaxIt; 
     
       if rand<0.5
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
            newPopPos(i,:)=PopPos(i,:)+S*(rand*BestX-rand*PopPos(i,:)); %Equation (8)
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
 
  
     HisBestFit(It)=BestF;
end


