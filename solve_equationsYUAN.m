%% ������⵱ǰ����������һʱ�����ݵĺ���
function x_true_final=solve_equations(var,target,correct_num,index_sum,index_num)
x_true_final=0; %����Ϊ�ս���������0�򷵻�0ֵ
%% ��������������һʱ�̵����ݵ��뷽�̵�������
P_i=1:index_sum; %��������������һʱ�̵����ݵ��뷽�̵�������
[varm,varn]=size(var);
re=1;
for ii =1:varn  
    re=re.*(var{ii}(correct_num).^(index_num(P_i,ii)));
end
P(P_i,1)=re;          
% P(P_i,1)=(b(correct_num).^(index_num(P_i,1))).*(c(correct_num).^(index_num(P_i,2))).*(d(correct_num).^(index_num(P_i,3)));
%% �õ�����������һʱ�̵�����
x_true=target*P; %���̵����������������ϵ����˵õ�����������һʱ�̵�����
%% �Է��̵Ľ���м���
if isempty(x_true)==1||isreal(x_true)==0||x_true==0 %�жϽ��Ƿ�Ϊ��,�Ƿ�Ϊʵ��,�Ƿ�Ϊ0  
return; %������һ���ж�Ϊ��,���������
end
x_true_final=x_true; %�õ���������Ľ�
end