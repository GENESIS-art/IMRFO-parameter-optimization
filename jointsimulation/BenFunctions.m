
function Fit=BenFunctions(X,FunIndex,Dim)

 global k1 k2 k3;

 switch FunIndex
     
   %%%%%%%%%%%%%%%%%%%%%%%%%%unimodal function%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  
%Sphere

%Rosenbrock
case 5

    load_system('pso');  % 加载模型

    param1 = X(1);  % 微分增益
    param2 = X(2);  % 比例增益
    param3 = X(3);  % 积分增益

    disp('Starting simulation...');
    disp(['Current parameters: Kp=', num2str(param2), ', Kd=', num2str(param1), ', Ki=', num2str(param3)]);
    fprintf('a\n');

    k1 = param1; 
    k2 = param2;
    k3 = param3;
    omega1 = 1;  % 根据需要设置的权重系数
    omega2 = 10;
    simOut = sim('IMRFO', 'ReturnWorkspaceOutputs', 'on', 'timeout', 100);
    error1 = simOut.error1;  % 获取误差信号
%     errorData = error1;  % 提取误差数据
%      % 计算误差
%     error = sum(errorData.^2);  % 使用平方和作为误差
%     Fit=error;
    e = error1;  % 误差信号
    t = (1:length(e)) * 0.001;  % 时间向量，Ts是采样时间间隔
    t = t';
    a = t .* abs(e);
    % 使用trapz进行数值积分，计算 |e(t)| 和 t|e(t)|
    integral1 = trapz(t, abs(e));  % |e(t)| 的积分
    integral2 = trapz(t, a);  % t|e(t)| 的积分
    
    % 误差函数
    J = omega1 * integral1 + omega2 * integral2;
    Fit=J;
    fprintf('Current PD Params: Kp = %f, Kd = %f, Ki = %f,Current Error = %f\n', k2, k1,k3, J);
    disp(J);
    fprintf('b\n');
%     close all;
    bdclose('pso');  % 每次仿真后关闭模型


end
 
