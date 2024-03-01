%% 用于求解当前滑动窗口下一时刻数据的函数
function x_true_final=solve_equations(var,target,correct_num,index_sum,index_num)
x_true_final=0; %若解为空解或虚数解或0则返回0值
%% 将关联变量在下一时刻的数据导入方程的所有项
P_i=1:index_sum; %将关联变量在下一时刻的数据导入方程的所有项
[varm,varn]=size(var);
re=1;
for ii =1:varn  
    re=re.*(var{ii}(correct_num).^(index_num(P_i,ii)));
end
P(P_i,1)=re;          
% P(P_i,1)=(b(correct_num).^(index_num(P_i,1))).*(c(correct_num).^(index_num(P_i,2))).*(d(correct_num).^(index_num(P_i,3)));
%% 得到待检测变量下一时刻的数据
x_true=target*P; %方程的所有项与所有项的系数相乘得到待检测变量下一时刻的数据
%% 对方程的解进行检验
if isempty(x_true)==1||isreal(x_true)==0||x_true==0 %判断解是否为空,是否为实数,是否为0  
return; %若其中一项判断为真,则结束程序
end
x_true_final=x_true; %得到经过检验的解
end