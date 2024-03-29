%% 稀疏自适应匹配追踪算法
function [theta]=SAMP(y,A,S)
[y_rows,y_columns]=size(y);
if y_rows<y_columns
y=y'; %y should be a column vector
end
[M,N]=size(A); %传感矩阵A为M*N矩阵
theta=zeros(N,1); %用来存储恢复的theta(列向量)
Pos_theta=[]; %用来迭代过程中存储A被选择的列序号
r_n=y; %初始化残差(residual)为y
L=S; %初始化步长(Size of the finalist in the first stage)
Stage=1; %初始化Stage
IterMax=M;
for i=1:IterMax %最多迭代M次
%% (1)Preliminary Test
product=A'*r_n; %传感矩阵A各列与残差的内积
[~,pos]=sort(abs(product),'descend'); %降序排列
Sk=pos(1:L); %选出最大的L个
%% (2)Make Candidate List
Ck=union(Pos_theta,Sk);
%% (3)Final Test
if length(Ck)<=M
At = A(:,Ck); %将A的这几列组成矩阵At
else
theta_ls=0;
break;
end
%% y=At*theta，以下求theta的最小二乘解(Least Square)
theta_ls=(At'*At)^(-1)*At'*y; %最小二乘解
[~,pos]=sort(abs(theta_ls),'descend'); %降序排列
F=Ck(pos(1:L));
%% (4)Compute Residue
%% A(:,F)*theta_ls是y在A(:,F)列空间上的正交投影
theta_ls=(A(:,F)'*A(:,F))^(-1)*A(:,F)'*y;
r_new=y-A(:,F)*theta_ls; %更新残差r
if norm(r_new)<1e-10 %halting condition true 
Pos_theta=F;
%% r_n=r_new，更新r_n以便输出最新的r_n
break; %quit the iteration
elseif norm(r_new)>=norm(r_n) %stage switching 
Stage=Stage+1; %Update the stage index 
L=Stage*S; %Update the size of finalist
if i==IterMax %最后一次循环
Pos_theta=F; %更新Pos_theta以与theta_ls匹配，防止报错
end
% i=i-1; %迭代次数不更新
else
Pos_theta=F; %Update the finalist Fk
r_n=r_new; %Update the residue
end
end
theta(Pos_theta)=theta_ls; %恢复出的theta
end