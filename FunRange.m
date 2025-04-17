function [Low,Up,Dim]=FunRange(FunIndex)
 
Dim=30;

switch FunIndex
    


    case 5
    Low=[0.5 1 0.5];Up=[300 5000 2000];Dim=3;




    otherwise 
    Low=0;Up=10;Dim=4;
    
end

