%% ϡ������Ӧƥ��׷���㷨
function [theta]=SAMP(y,A,S)
[y_rows,y_columns]=size(y);
if y_rows<y_columns
y=y'; %y should be a column vector
end
[M,N]=size(A); %���о���AΪM*N����
theta=zeros(N,1); %�����洢�ָ���theta(������)
Pos_theta=[]; %�������������д洢A��ѡ��������
r_n=y; %��ʼ���в�(residual)Ϊy
L=S; %��ʼ������(Size of the finalist in the first stage)
Stage=1; %��ʼ��Stage
IterMax=M;
for i=1:IterMax %������M��
%% (1)Preliminary Test
product=A'*r_n; %���о���A������в���ڻ�
[~,pos]=sort(abs(product),'descend'); %��������
Sk=pos(1:L); %ѡ������L��
%% (2)Make Candidate List
Ck=union(Pos_theta,Sk);
%% (3)Final Test
if length(Ck)<=M
At = A(:,Ck); %��A���⼸����ɾ���At
else
theta_ls=0;
break;
end
%% y=At*theta��������theta����С���˽�(Least Square)
theta_ls=(At'*At)^(-1)*At'*y; %��С���˽�
[~,pos]=sort(abs(theta_ls),'descend'); %��������
F=Ck(pos(1:L));
%% (4)Compute Residue
%% A(:,F)*theta_ls��y��A(:,F)�пռ��ϵ�����ͶӰ
theta_ls=(A(:,F)'*A(:,F))^(-1)*A(:,F)'*y;
r_new=y-A(:,F)*theta_ls; %���²в�r
if norm(r_new)<1e-10 %halting condition true 
Pos_theta=F;
%% r_n=r_new������r_n�Ա�������µ�r_n
break; %quit the iteration
elseif norm(r_new)>=norm(r_n) %stage switching 
Stage=Stage+1; %Update the stage index 
L=Stage*S; %Update the size of finalist
if i==IterMax %���һ��ѭ��
Pos_theta=F; %����Pos_theta����theta_lsƥ�䣬��ֹ����
end
% i=i-1; %��������������
else
Pos_theta=F; %Update the finalist Fk
r_n=r_new; %Update the residue
end
end
theta(Pos_theta)=theta_ls; %�ָ�����theta
end