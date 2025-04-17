 global k1 k2 k3;
 MaxIteration=20; 
 PopSize=8;
 FunIndex=5;

[BestX,BestF,HisBestF]=IMRFO(FunIndex,MaxIteration,PopSize);



display(['F_index=', num2str(FunIndex)]);
display(['The best fitness is: ', num2str(BestF)]);
Optimal(FunIndex)=BestF;

%display(['The best solution is: ', num2str(BestX)]);
 if BestF>=0
     semilogy(HisBestF,'r','LineWidth',2);
 else
     plot(HisBestF,'r','LineWidth',2);
 end
 xlabel('Iterations');
 ylabel('Fitness');
 title(['F',num2str(FunIndex)]);








