%% ����ѹ����֪�����̹�ҵ�쳣�������ʵʱ�������޸��㷨���޻���ѵ�����̣�
%% �����ر�
clear;close all;clc; %��չ��������ر����д��ڣ������������
%% ��ʼ��ʱ
tic; %��ʼ��ʱ
%% ��������
[Y1,Y2,Y3,Y4]=coupled_logistic_4(0.4,0.4,0.4,0.4,10000); %���������ļ�
a(:,1)=Y1(2000:5000); %�����������a������
b(:,1)=Y1(1999:4999); %�����������b������
c(:,1)=Y3(1999:4999); %�����������c������
% d(:,1)=Y2; %�����������c������
var={a,b,c};
var_bcd={b,c};
[varm,varn]=size(var);
%% ��֪�����ָ������
index_num = Combination_MN([3,3]);
[index_sum,colum]=size(index_num); %����ָ������������ϵ������
% index_num_i=0; %����ָ����������ʼֵ
% index_num=zeros(index_sum,2); %��������ָ�����
% for index_i=0:1:2 %ѭ���õ����е�ָ�����
%     for index_j=0:1:2
% %         for index_k=0:1:2
%             index_num_i=index_num_i+1;
%             index_num(index_num_i,:)=[index_i,index_j]; %�õ����е�ָ�����
% %         end
%     end
% end
%% �����������a��ʵ������
aa=a; %�����������a��ʵ������
%% �ɵ���������
rebuild_num=1000; %���ñ���a�Ĵ�������������
%% ����(1)
slip_start=2000; %���ô��ڻ�������ʼ��
slip_end=slip_start+rebuild_num-1; %���ô��ڻ����Ľ�ֹ��
%% ����(2)
coefficient_num=index_sum; %����ϵ������
%% ����(3)������ǿ��أ�
judge_num=5; %�����жϴ������ƣ��м��ѡȡλ���ɿ�ѡ�޸����������(ԭΪ10)
alternative_relative_error_judge_limit=1e-2; %���ÿ�ѡ�޸����������ƣ�ԭΪ1e-2��
%% ����(4)������ǿ��أ����ع��ɹ��������£�ÿ��ѡ������������ܵõ�һ������ƽ���ع����������Ƶ�һ�飩
loop_num_limit=5; %���ô�ѡ����������Ŀ�����������(ԭΪ10)
choose_num=5; %����ϸѡ�����������趨����������(ԭΪ10)
loop_num_coefficient=2; %���ò�������������
%% ����(5)�������޸������������Ӧǿ��ƽ���ع�����������ƣ�
mean_error_limit=1e-2; %����ƽ���ع����������ƣ�ԭΪ1e-2��
most_mean_error_limit=1e-3; %���ô󲿷��ع�����������
%% ����(6)�������޸������������Ӧǿ��ƽ���ع�����������ƣ�
%% �ɵ�
test_error_judge_limit=1e-3; %���ò��Ե���޸����������ƣ�ԭΪ1e-3��
%% ����(7)�������޸������������Ӧǿ��ƽ���ع�����������ƣ�
middle_relative_error_judge_limit=1e-3; %�������÷���ʱ���޸����������ƣ�ԭΪ1e-3��
%% ����(8)�������޸������������Ӧǿ��ƽ���ع�����������ƣ�
next_relative_error=1e-3; %������һ���޸����������ƣ�ԭΪ1e-3��
error_expansion_coefficient=1; %������һ���޸�����������ϵ��(ԭΪ1������ϵ����Ӧ�ù����ô����ʵ����Ч�ʣ�
%% ����(10)
sampling_num=8; %���ò�����ĸ�����ʼֵΪ1���������0��
sampling_range=50; %���ò�����Χ��ʼֵΪ10���������0��
%% ����ȫ������
final_target=zeros(rebuild_num,coefficient_num); %�������е�ϵ������
all_loop_num=zeros(rebuild_num,1); %�������е�ѭ������
all_judge_times=zeros(rebuild_num,1); %�������е��жϴ���
all_alternative_a_correct=zeros(judge_num,rebuild_num); %�������еı�ѡ�޸�����
all_alternative_relative_error_judge=zeros(judge_num,rebuild_num); % �������еı�ѡ�޸�������
all_test_error=zeros(choose_num,rebuild_num); % �������еĲ��Ե��޸�������
all_only_minposition=zeros(rebuild_num,1); %�������еĵ�ǰ��С�޸��������λ�ã�����(4)��
all_sampling_num=zeros(rebuild_num,1); %�������еĲ��������
all_sampling_range=zeros(rebuild_num,1); %�������еĲ�����Χ
%% �����һ�γ��ֵķ��̼�λ�ã���final_target�У�
first_equations_combination=zeros(rebuild_num,coefficient_num); %�����һ�γ��ֵķ���
first_equations_combination_local=zeros(rebuild_num,1); %�����һ�γ��ֵķ���λ��
%% ���淽��ϵ������
target=zeros(1,coefficient_num); %���淽��ϵ������
%% ����������
hwait=waitbar(0,'please wait...'); %���������
set(hwait,'doublebuffer','on'); %����˫���棬��������̬��ʾ�ٷֱ�ʱ����������
step=(rebuild_num)/100; %�������ܳ���Ϊ100ʱ��Ӵ���
%% ����ɹ�����һʱ�̴��ݵ��м�ϵ��������ѭ�������Ч�ʺ��ȶ��ԣ�
record_local=zeros(rebuild_num,1); %��¼����֮ǰ���̵����ݵ�λ��
record_error=zeros(rebuild_num,1); %��¼����֮ǰ���̵��޸��������С
%% �����������ζԼ�����ݽ��м������޸�
for slip=slip_start:1:slip_end %���ڻ�����ʼ
%% ����ÿ�λ�������Ҫ���õĲ���
    relative_error_judge_limit=next_relative_error; %������һ���޸�����������
%% ���汸ѡ����
    alternative_a_correct=zeros(judge_num,1); %���汸ѡ���޸�����
    alternative_relative_error_judge=zeros(judge_num,1); %���汸ѡ���޸�������
    alternative_final_target=zeros(judge_num,coefficient_num); %���汸ѡ�ķ���ϵ������
%% ������ʼ�жϴ���
    judge_times=0; %������ʼ�жϴ���
%% ����������
    PerStr=fix((slip-slip_start+1)/step);
    str=['����������',num2str(PerStr),'%'];
    waitbar((slip-slip_start+1)/(rebuild_num),hwait,str); %��ɽ���������
%% �õ�����Ҫ��ķ���,����ѭ��������֮ǰ����
    test_num=slip; %���ò��Ե�λ��
    for f=(slip-slip_start):-1:1 %���γ���֮ǰ���е�һ�γ��ֵķ���
        if first_equations_combination_local(f,1)~=0 %����÷����ǵ�һ�γ��֣����������ж�
%% �õ���һʱ�̵��޸�����
            target=final_target(first_equations_combination_local(f,1),:); %�õ�����ϵ��
            a_correct=solve_equationsYUAN(var_bcd,target,test_num+1,index_sum,index_num); %�õ���һʱ�̵��޸�����
%% �ж��޸������Ƿ����Ҫ����������Ҫ���������޸����ݴ����������Ա�ɸѡ
                if a_correct~=0 %��Ҫ�жϣ�����һʱ�̵��޸����ݲ�Ϊ0���������������
                    relative_error_judge=abs((a_correct-a(test_num+1))/a_correct); %��һʱ�̵��޸������������һʱ�̵������޷�ȷ���Ƿ��쳣�����Է�ĸΪ�޸����ݣ�
                    record_error(1+slip-slip_start,1)=relative_error_judge; %��¼����֮ǰ���̵��޸�������
                        if relative_error_judge<middle_relative_error_judge_limit %���޸�������С����ֵ��˵���޸���������Ҫ�󣬼����ж���һʱ�̵��޸�����
                            final_target(1+slip-slip_start,:)=target; %�������յ�ϵ������
                            a(slip+1)=a_correct; %ֱ�ӽ��޸����ݴ��룬����޸�
% a(slip+1)=a(test_num+1); %�������쳣�������޸�
                            record_local(slip-slip_start+1,1)=first_equations_combination_local(f,1); %��¼���õ�����һ���ķ��̣�slip-slip_start+1����ֹ��0������������
                            break; %������ǰѭ��
                        end
                end
        end
    end 
    if record_local(slip-slip_start+1,1)~=0 %�����ǰ������֮ǰ�ķ����򴰿ڼ�������
        continue; %������һѭ��
    end
%% ����ȫ������
    all_sampling_num(slip-slip_start+1,1)=sampling_num; %�������еĲ��������
    all_sampling_range(slip-slip_start+1,1)=sampling_range; %�������еĲ�����Χ
%% ������ʼѭ������
    loop_num_start=loop_num_limit*loop_num_coefficient; %������ʼѭ������
%% ������Ե��޸��������ж�ѭ��
    test_error_judge=1; %���ò��Ե��޸��������ʼֵ������ѭ��
    while(test_error_judge>=test_error_judge_limit) %�����Ե��޸�������С����ֵ������ѭ��
%% ����ƽ���ع����������Ҫ��ķ��������ж�ѭ��
        all_mean_error_pass_num=0; %����ƽ���ع����������Ҫ��ķ���������ʼֵ������ѭ��
        loop_num=loop_num_start; %������ʼѭ������
        while(all_mean_error_pass_num<loop_num_limit) %��ƽ���ع����������Ҫ��ķ����������ڵ�����ֵ������ѭ��
%% ��̬����ѭ������
            loop_num=loop_num+1; %�������������㣬ѭ������+1
%% ����ÿ��ѭ���õ�������
            all_target=zeros(loop_num,coefficient_num); %����ÿ��ѭ��������ϵ������
            all_sample=zeros(loop_num,sampling_num); %����ÿ��ѭ�������в�����
            all_mean_error=zeros(loop_num,1); %����ÿ��ѭ��������ƽ���ع�������
            all_a_error=zeros(loop_num,1);
%% ѭ����ȡƽ���ع�������
            for i=1:loop_num %�涨����ѭ����ʼ
%% �������������
                sample=randperm(sampling_range,sampling_num); %������Χ������ɼ��涨�����Ĳ�����
                sample=sort(sample); %�������С��������
%% �õ���������
                Phi=zeros(sampling_num,coefficient_num); %�����������
                Phi_row=zeros(1,index_sum); %������������������
                for Phi_i=1:sampling_num %ѭ���õ��������������������
%% �õ���������������������Ż�forѭ����
                    Phi_row_i=1:index_sum;
                    re=1;
                    for ii =2:varn  
                        re=re.*(var{ii}(sample(Phi_i)+slip-sampling_range).^(index_num(Phi_row_i,ii-1)));
                    end
                    Phi_row(1,Phi_row_i)=re;
%                     Phi_row(1,Phi_row_i)=(b(sample(Phi_i)+slip-sampling_range).^(index_num(Phi_row_i,1))).*(c(sample(Phi_i)+slip-sampling_range).^(index_num(Phi_row_i,2))).*(d(sample(Phi_i)+slip-sampling_range).^(index_num(Phi_row_i,3)));
%% �õ���������
                    Phi(Phi_i,:)=Phi_row; %�õ���������
                end
%% �õ�Phi�������������L2���������Ż�forѭ����
                l=zeros(coefficient_num,1); %����Phi�������������L2����
                s=1:coefficient_num;
                l(s)=(sum(Phi(:,s).^2)).^0.5; %�õ�Phi�������������L2����
%% ���������й�һ�� 
                A=Phi./repmat(sqrt(sum(Phi.^2,1)),size(Phi,1),1); %���������й�һ��
%% �õ��۲�����y�����Ż�forѭ����
                y=zeros(sampling_num,1); %����۲����y
                y_i=1:sampling_num;
                y(y_i)=a(sample(y_i)+slip-sampling_range); %�õ��۲�����y
%% ����SAMP�����õ�����ϵ��
                theta=SAMP(y,A,1); %����SAMP�����õ�����ϵ��
%% ���һ�����õ�����ϵ�������Ż�forѭ����
                coefficient=zeros(1,coefficient_num); %���淽��ϵ��
                target_samp_i=1:coefficient_num;
                coefficient(target_samp_i)=theta(target_samp_i)./l(target_samp_i); %���һ�����õ�����ϵ��
%% ����ÿ��ѭ��������ϵ�����󼰲�����
                all_target(i,:)=coefficient; %����ÿ��ѭ��������ϵ������
                all_sample(i,:)=sample; %����ÿ��ѭ�������в�����
%% �õ��ع����ݣ����Ż�forѭ����
                a_rebuild_resolution=zeros(sampling_range,index_sum); %�����������ϵ����Ĳ�������
                h=1:sampling_range; %ѭ���õ����������ڵ��ع�����
%% �õ���������ϵ����Ĳ�������
                for a_rebuild_i=1:index_sum
                    re2=1;
                    for jjj =2:varn  
                        re2=re2.*(var{jjj}(h+slip-sampling_range).^(index_num(a_rebuild_i,jjj-1)));
                    end
                    a_rebuild_resolution(:,a_rebuild_i)=coefficient(1,a_rebuild_i)*re2;
%                     a_rebuild_resolution(:,a_rebuild_i)=coefficient(1,a_rebuild_i)*(b(h+slip-sampling_range).^(index_num(a_rebuild_i,1))).*(c(h+slip-sampling_range).^(index_num(a_rebuild_i,2))).*(d(h+slip-sampling_range).^(index_num(a_rebuild_i,3)));
                end
                a_rebuild=sum(a_rebuild_resolution,2); %�����ع�����
%% �õ�ƽ���ع�������
                a_real=a((slip-sampling_range+1):(slip),1); %�õ�������Χ�ڸ������ʵ������
                a_error=(a_rebuild-a_real)./a_real; %�õ�������Χ�ڸ�������ع���������ڴ����ڵ�����ȷ�������쳣�ģ����Է�ĸΪʵ�����ݣ�
                all_mean_error(i,1)=sum(abs(a_error))/sampling_range; %��������ƽ���ع������ƽ���ع���������������ع����������ֵ֮�͵�ƽ��ֵ
                if sum(abs(a_error)<mean_error_limit)==(sampling_range)
                    all_a_error(i,1)=1;
                end
            end
%% �ҵ�����ƽ���ع����������Ҫ��ķ���λ�ü�����
            all_mean_error_pass=find(all_a_error==1); %�ҵ�����ƽ���ع����������Ҫ��ķ���λ��
            all_mean_error_pass_num=length(all_mean_error_pass); %�ҵ�����ƽ���ع����������Ҫ��ķ��̸���
        end
%% ѡ��ƽ���ع���������С�����ɸ�����
        [result,index]=sort(all_mean_error); %ƽ���ع�����������
        all_mean_error_pass_part=index(1:choose_num,1); %ѡ�񷽳̲ɼ�����������ƽ���ع���������С����ѡ��
%% �õ����о���ѡ��ķ����ڲ��Ե���޸������Լ��޸�������
        test_num=slip; %���ò��Ե�λ��
        test_error=zeros(choose_num,1); %������Ե��޸�������
        for r=1:choose_num %ѭ��������о���ѡ��ķ����ڲ��Ե���޸�����
        target=all_target(all_mean_error_pass_part(r,1),:); %�õ�����ϵ��
        test_rebuild_01=solve_equationsYUAN(var_bcd,target,test_num,index_sum,index_num); %�õ����Ե���޸�����
        test_real_01=a(test_num,1); %�õ��������ݵ�ʵ������
        test_error_01=abs((test_rebuild_01-test_real_01)/test_real_01); %������Ե���޸���������ڲ��Ե������ȷ�������쳣�ģ����Է�ĸΪʵ�����ݣ�
        test_error(r,1)=test_error_01; %������Ե���޸���������������ݵ��޸�������ʵ�����ݵ�������ľ���ֵ��
        end
%% �������еĲ��Ե��޸�������
        all_test_error(:,slip-slip_start+1)=test_error; %�������еĲ��Ե��޸�������
%% ������while�����µõ���һʱ�̵��޸�����
        a_correct=0; %����a_correct��ʼֵ������ѭ��
        test_error_judge=0; %����test_error_judge��ʼֵ������ѭ��
        while(a_correct==0&&test_error_judge<test_error_judge_limit) %����һʱ�̵��޸����ݲ�Ϊ�����С�������еĲ��Ե��޸���������ڵ�����ֵ��������ѭ��
            minposition=find(test_error==min(test_error)); %�ҵ���ǰ��С�Ĳ��Ե��޸��������λ��
            only_minposition=minposition(1,1); %�õ���ǰΨһ����С�޸��������λ��
            test_error_judge=test_error(only_minposition); %�õ���ǰ��С�Ĳ��Ե��޸�������
            test_error(minposition)=1; %��ǰ����С���Ե��޸���������1
            if test_error_judge<test_error_judge_limit  %��Ҫ�жϣ�����ǰ��С�޸�������С����ֵ���޸���һʱ�̵�����
                target=all_target(all_mean_error_pass_part(only_minposition,1),:); %�õ�����ϵ��
                a_correct=solve_equationsYUAN(var_bcd,target,test_num+1,index_sum,index_num); %�õ���һʱ�̵��޸�����
            end
        end
%% �������еĵ�ǰ��С�޸��������λ��
        all_only_minposition(slip-slip_start+1,1)=only_minposition; %�������еĵ�ǰ��С�޸��������λ��
%% ��ֹѭ���������
        loop_num=loop_num-1; %��һ�ν����ж�ƽ���ع����������Ҫ��ķ��������Ƿ��㹻��ѭ��ʱ����֤ѭ���������䣨��һ��ѭ���Ѿ�֤����ǰѭ�������㹻��
%% �ж��޸������Ƿ����Ҫ����������Ҫ���������޸����ݴ����������Ա�ɸѡ
%% ��Ҫ�жϣ�����һʱ�̵��޸����ݲ�Ϊ0���������������
        if a_correct~=0
            relative_error_judge=abs((a_correct-a(test_num+1))/a_correct); %��һʱ�̵��޸������������һ������ݲ�ȷ���Ƿ��쳣�ģ����Է�ĸΪ�޸����ݣ�
%% �������������ֵ�����жϴ���С����ֵ��������������̣����жϴ�������ֵ�ڵõ�С����ֵ���޸������˵���޸���������Ҫ��
            if relative_error_judge>=relative_error_judge_limit&&judge_times<judge_num 
%% ���汸ѡ���޸������Լ��޸�������
                alternative_a_correct(judge_times+1,1)=a_correct; %���汸ѡ���޸�����
                alternative_relative_error_judge(judge_times+1,1)=relative_error_judge; %���汸ѡ���޸�������
%% ���汸ѡ��ϵ�����󼰲�����
                alternative_final_target(judge_times+1,:)=target; %���汸ѡ��ϵ������
%% ����ѭ�����ع�����ѭ��
                judge_times=judge_times+1; %�жϴ���+1
                relative_error_judge_limit=relative_error_judge_limit*error_expansion_coefficient; %��С��һ���޸����������ƣ�*�޸����������ϵ����
                test_error_judge=1; %Ϊcontinue���½���ѭ�����ó�ֵ
                continue; %���½����жϲ��Ե��޸��������Ƿ����Ҫ���ѭ����
            end
        end
    end
%% �����޸����̣�1��
    if judge_times<judge_num %���жϴ���С����ֵ��˵���ڹ涨�����ڵõ��˷���Ҫ����޸�����
        a(slip+1)=a_correct; %ֱ�ӽ��޸����ݴ��룬����޸�
% a(slip+1)=a(test_num+1); %�������쳣�������޸�
%% �������յ�ϵ������
        final_target(1+slip-slip_start,:)=target; %�������յ�ϵ������
%% �����һ�γ��ֵķ��̼�λ�ã���final_target�У�
        first_equations_combination(1+slip-slip_start,:)=target; %�����һ�γ��ֵķ��̣���final_target�У�
        first_equations_combination_local(1+slip-slip_start,1)=1+slip-slip_start; %�����һ�γ��ֵķ���λ�ã���final_target�У�
    elseif judge_times==judge_num  %���жϴ���������ֵ��˵���ڹ涨������û�еõ�����Ҫ����޸����ݣ�����ɸѡ����
        [j_result,j_index]=sort(alternative_relative_error_judge); %���޸������������������
%% ѡȡ�м�㣨�����쳣�Ƿ�ȷ���������޸����ݴ��룬����޸���
        if sum(alternative_relative_error_judge>=alternative_relative_error_judge_limit)==judge_num
            middle_num=ceil(judge_num/2); %�����п�ѡ�޸���������ڻ������ֵ������������������˵��������ݺܿ��ܳ����쳣���м��ȡ�м�ֵ
        elseif sum(alternative_relative_error_judge>=alternative_relative_error_judge_limit)<judge_num
            middle_num=1; %�����п�ѡ�޸���������ڻ������ֵ������С����������˵��������ݿ������쳣�����������ϴ��ԭ������������������м��ȡ��Сֵ    
        end
%% �����޸����̣�2��
        find_min_num=j_index(middle_num,1); %ѡ���ŵڼ�С���޸�������õ���λ��
        a_correct=alternative_a_correct(find_min_num,1); %�õ���ѡ���޸�����
        a(slip+1)=a_correct; %���޸����ݴ��룬����޸�
%% �������յ�ϵ������
        final_target(1+slip-slip_start,:)=alternative_final_target(find_min_num,:); %�������յ�ϵ������
%% �����һ�γ��ֵķ��̼�λ�ã���final_target�У�
        first_equations_combination(1+slip-slip_start,:)=alternative_final_target(find_min_num,:); %�����һ�γ��ֵķ���
        first_equations_combination_local(1+slip-slip_start,1)=1+slip-slip_start; %�����һ�γ��ֵķ���λ�ã���final_target��
    end
%% ���浱ǰ��ѡ���޸�����
    all_alternative_a_correct(:,1+slip-slip_start)=alternative_a_correct; % ���浱ǰ��ѡ���޸�����
%% ���浱ǰ��ѡ���޸�������
    all_alternative_relative_error_judge(:,1+slip-slip_start)=alternative_relative_error_judge; % ���浱ǰ��ѡ���޸�������
%% ���浱ǰѭ���������жϴ���
    loop_num=loop_num+1; %�ָ�ѭ������
    all_loop_num(1+slip-slip_start,1)=loop_num; %���浱ǰѭ������
    all_judge_times(1+slip-slip_start,1)=judge_times; %���浱ǰ�жϴ���
end
%% �������ر�
close(hwait); %�������ر�%% ʵ�����ݺ��޸����ݴ���
%% ʵ�����ݺ��޸����ݴ���
x_final(1,:)=aa((slip_start+1):(slip_end+1)); %�õ�ʵ������ 
x_final(2,:)=a((slip_start+1):(slip_end+1)); %�õ��޸�����
normal_repair(1,:)=x_final(2,:)-x_final(1,:); %�õ�ʵ���������޸����ݵ����
relative_normal_repair(1,:)=normal_repair(1,:)./x_final(1,:); %�õ�ʵ���������޸����ݵ���������ڼ����㷨����֪ʵ�����ݶ������쳣�ģ����޸����ݾ�����Σ�
relative_normal_repair(2,:)=normal_repair(1,:)./x_final(2,:); %�õ�ʵ���������޸����ݵ�����������޸����ݸ��ɿ������Է�ĸΪ�޸����ݣ�
%% ʵ���������޸����ݵĶԱ�ͼ
figure; %��ͼ��ʼ
plot(1:rebuild_num,x_final(1,:),'ks-','MarkerSize',3,'MarkerFaceColor','k'); %����ʵ����������
hold on; %ͼ��ϲ�
plot(1:rebuild_num,x_final(2,:),'r^-','MarkerSize',3,'MarkerFaceColor','r'); %�����޸���������
% axis([1 rebuild_num 1.10 1.40]); %���������᷶Χ
legend('ʵ������','�޸�����'); %ͼ���ע
hold off; %��ͼ����
%% �޸�����������
figure; %��ͼ��ʼ
plot(1:rebuild_num,relative_normal_repair(1,:),'ko-','MarkerSize',3,'MarkerFaceColor','k'); %�����޸�����������
hold off; %��ͼ����
%% ������ʱ
toc; %������ʱ