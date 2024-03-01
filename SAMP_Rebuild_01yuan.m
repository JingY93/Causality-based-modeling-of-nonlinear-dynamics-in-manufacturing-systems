%% 采用压缩感知的流程工业异常监测数据实时检验与修复算法（无滑动训练过程）
%% 清空与关闭
clear;close all;clc; %清空工作区，关闭所有窗口，清空命令区域
%% 开始计时
tic; %开始计时
%% 导入数据
[Y1,Y2,Y3,Y4]=coupled_logistic_4(0.4,0.4,0.4,0.4,10000); %导入数据文件
a(:,1)=Y1(2000:5000); %导入待检测变量a的数据
b(:,1)=Y1(1999:4999); %导入关联变量b的数据
c(:,1)=Y3(1999:4999); %导入关联变量c的数据
% d(:,1)=Y2; %导入关联变量c的数据
var={a,b,c};
var_bcd={b,c};
[varm,varn]=size(var);
%% 感知矩阵的指数构成
index_num = Combination_MN([3,3]);
[index_sum,colum]=size(index_num); %设置指数的总组数即系数个数
% index_num_i=0; %设置指数组数的起始值
% index_num=zeros(index_sum,2); %储存所有指数组合
% for index_i=0:1:2 %循环得到所有的指数组合
%     for index_j=0:1:2
% %         for index_k=0:1:2
%             index_num_i=index_num_i+1;
%             index_num(index_num_i,:)=[index_i,index_j]; %得到所有的指数组合
% %         end
%     end
% end
%% 保存待检测变量a的实际数据
aa=a; %保存待检测变量a的实际数据
%% 可调参数设置
rebuild_num=1000; %设置变量a的待检验数据总量
%% 参数(1)
slip_start=2000; %设置窗口滑动的起始点
slip_end=slip_start+rebuild_num-1; %设置窗口滑动的截止点
%% 参数(2)
coefficient_num=index_sum; %设置系数个数
%% 参数(3)（精度强相关）
judge_num=5; %设置判断次数限制，中间点选取位置由可选修复相对误差决定(原为10)
alternative_relative_error_judge_limit=1e-2; %设置可选修复相对误差限制（原为1e-2）
%% 参数(4)（精度强相关，在重构成功率限制下，每次选择几组采样组数能得到一定满足平均重构相对误差限制的一组）
loop_num_limit=5; %设置粗选采样组数（目标采样组数）(原为10)
choose_num=5; %设置细选采样组数（设定采样组数）(原为10)
loop_num_coefficient=2; %设置采样组数扩大倍数
%% 参数(5)（单点修复相对误差的限制应强于平均重构相对误差的限制）
mean_error_limit=1e-2; %设置平均重构相对误差限制（原为1e-2）
most_mean_error_limit=1e-3; %设置大部分重构相对误差限制
%% 参数(6)（单点修复相对误差的限制应强于平均重构相对误差的限制）
%% 可调
test_error_judge_limit=1e-3; %设置测试点的修复相对误差限制（原为1e-3）
%% 参数(7)（单点修复相对误差的限制应强于平均重构相对误差的限制）
middle_relative_error_judge_limit=1e-3; %设置沿用方程时的修复相对误差限制（原为1e-3）
%% 参数(8)（单点修复相对误差的限制应强于平均重构相对误差的限制）
next_relative_error=1e-3; %设置下一点修复相对误差限制（原为1e-3）
error_expansion_coefficient=1; %设置下一点修复相对误差扩张系数(原为1，扩张系数不应该过大，用处是适当提高效率）
%% 参数(10)
sampling_num=8; %设置采样点的个数初始值为1（必须大于0）
sampling_range=50; %设置采样范围初始值为10（必须大于0）
%% 储存全局数据
final_target=zeros(rebuild_num,coefficient_num); %储存所有的系数矩阵
all_loop_num=zeros(rebuild_num,1); %储存所有的循环次数
all_judge_times=zeros(rebuild_num,1); %储存所有的判断次数
all_alternative_a_correct=zeros(judge_num,rebuild_num); %储存所有的备选修复数据
all_alternative_relative_error_judge=zeros(judge_num,rebuild_num); % 储存所有的备选修复相对误差
all_test_error=zeros(choose_num,rebuild_num); % 储存所有的测试点修复相对误差
all_only_minposition=zeros(rebuild_num,1); %储存所有的当前最小修复相对误差的位置（参数(4)）
all_sampling_num=zeros(rebuild_num,1); %储存所有的采样点个数
all_sampling_range=zeros(rebuild_num,1); %储存所有的采样范围
%% 储存第一次出现的方程及位置（在final_target中）
first_equations_combination=zeros(rebuild_num,coefficient_num); %储存第一次出现的方程
first_equations_combination_local=zeros(rebuild_num,1); %储存第一次出现的方程位置
%% 储存方程系数矩阵
target=zeros(1,coefficient_num); %储存方程系数矩阵
%% 进度条设置
hwait=waitbar(0,'please wait...'); %引入进度条
set(hwait,'doublebuffer','on'); %启用双缓存，以消除动态显示百分比时的闪屏问题
step=(rebuild_num)/100; %进度条总长不为100时添加此项
%% 储存可供向下一时刻传递的中间系数（减少循环，提高效率和稳定性）
record_local=zeros(rebuild_num,1); %记录沿用之前方程的数据点位置
record_error=zeros(rebuild_num,1); %记录沿用之前方程的修复相对误差大小
%% 滑动窗口依次对监测数据进行检验与修复
for slip=slip_start:1:slip_end %窗口滑动开始
%% 设置每次滑动后都需要重置的参数
    relative_error_judge_limit=next_relative_error; %设置下一点修复相对误差限制
%% 储存备选数据
    alternative_a_correct=zeros(judge_num,1); %储存备选的修复数据
    alternative_relative_error_judge=zeros(judge_num,1); %储存备选的修复相对误差
    alternative_final_target=zeros(judge_num,coefficient_num); %储存备选的方程系数矩阵
%% 设置起始判断次数
    judge_times=0; %设置起始判断次数
%% 进度条调用
    PerStr=fix((slip-slip_start+1)/step);
    str=['正在运行中',num2str(PerStr),'%'];
    waitbar((slip-slip_start+1)/(rebuild_num),hwait,str); %完成进度条调用
%% 得到符合要求的方程,不用循环，沿用之前方程
    test_num=slip; %设置测试点位置
    for f=(slip-slip_start):-1:1 %依次尝试之前所有第一次出现的方程
        if first_equations_combination_local(f,1)~=0 %如果该方程是第一次出现，进行如下判断
%% 得到下一时刻的修复数据
            target=final_target(first_equations_combination_local(f,1),:); %得到方程系数
            a_correct=solve_equationsYUAN(var_bcd,target,test_num+1,index_sum,index_num); %得到下一时刻的修复数据
%% 判断修复数据是否符合要求，若不符合要求则将所有修复数据储存起来，以备筛选
                if a_correct~=0 %重要判断，若下一时刻的修复数据不为0，则进入以下流程
                    relative_error_judge=abs((a_correct-a(test_num+1))/a_correct); %下一时刻的修复相对误差（由于下一时刻的数据无法确认是否异常，所以分母为修复数据）
                    record_error(1+slip-slip_start,1)=relative_error_judge; %记录沿用之前方程的修复相对误差
                        if relative_error_judge<middle_relative_error_judge_limit %若修复相对误差小于阈值，说明修复数据满足要求，继续判断下一时刻的修复数据
                            final_target(1+slip-slip_start,:)=target; %储存最终的系数矩阵
                            a(slip+1)=a_correct; %直接将修复数据代入，完成修复
% a(slip+1)=a(test_num+1); %数据无异常，不用修复
                            record_local(slip-slip_start+1,1)=first_equations_combination_local(f,1); %记录延用的是哪一步的方程（slip-slip_start+1，防止第0步看不出来）
                            break; %跳出当前循环
                        end
                end
        end
    end 
    if record_local(slip-slip_start+1,1)~=0 %如果当前延用了之前的方程则窗口继续滑动
        continue; %进入下一循环
    end
%% 储存全局数据
    all_sampling_num(slip-slip_start+1,1)=sampling_num; %储存所有的采样点个数
    all_sampling_range(slip-slip_start+1,1)=sampling_range; %储存所有的采样范围
%% 设置起始循环次数
    loop_num_start=loop_num_limit*loop_num_coefficient; %设置起始循环次数
%% 进入测试点修复相对误差判断循环
    test_error_judge=1; %设置测试点修复相对误差初始值，进入循环
    while(test_error_judge>=test_error_judge_limit) %若测试点修复相对误差小于阈值则跳出循环
%% 进入平均重构相对误差符合要求的方程数量判断循环
        all_mean_error_pass_num=0; %设置平均重构相对误差符合要求的方程数量初始值，进入循环
        loop_num=loop_num_start; %设置起始循环次数
        while(all_mean_error_pass_num<loop_num_limit) %若平均重构相对误差符合要求的方程数量大于等于阈值则跳出循环
%% 动态调整循环次数
            loop_num=loop_num+1; %若方程数量不足，循环次数+1
%% 储存每次循环得到的数据
            all_target=zeros(loop_num,coefficient_num); %储存每次循环的所有系数矩阵
            all_sample=zeros(loop_num,sampling_num); %储存每次循环的所有采样点
            all_mean_error=zeros(loop_num,1); %储存每次循环的所有平均重构相对误差
            all_a_error=zeros(loop_num,1);
%% 循环求取平均重构相对误差
            for i=1:loop_num %规定次数循环开始
%% 生成随机采样点
                sample=randperm(sampling_range,sampling_num); %采样范围内随机采集规定个数的采样点
                sample=sort(sample); %采样点从小到大排序
%% 得到测量矩阵
                Phi=zeros(sampling_num,coefficient_num); %储存测量矩阵
                Phi_row=zeros(1,index_sum); %储存测量矩阵的行向量
                for Phi_i=1:sampling_num %循环得到测量矩阵的所有行向量
%% 得到测量矩阵的行向量（已优化for循环）
                    Phi_row_i=1:index_sum;
                    re=1;
                    for ii =2:varn  
                        re=re.*(var{ii}(sample(Phi_i)+slip-sampling_range).^(index_num(Phi_row_i,ii-1)));
                    end
                    Phi_row(1,Phi_row_i)=re;
%                     Phi_row(1,Phi_row_i)=(b(sample(Phi_i)+slip-sampling_range).^(index_num(Phi_row_i,1))).*(c(sample(Phi_i)+slip-sampling_range).^(index_num(Phi_row_i,2))).*(d(sample(Phi_i)+slip-sampling_range).^(index_num(Phi_row_i,3)));
%% 得到测量矩阵
                    Phi(Phi_i,:)=Phi_row; %得到测量矩阵
                end
%% 得到Phi矩阵各列向量的L2范数（已优化for循环）
                l=zeros(coefficient_num,1); %储存Phi矩阵各列向量的L2范数
                s=1:coefficient_num;
                l(s)=(sum(Phi(:,s).^2)).^0.5; %得到Phi矩阵各列向量的L2范数
%% 测量矩阵按列归一化 
                A=Phi./repmat(sqrt(sum(Phi.^2,1)),size(Phi,1),1); %测量矩阵按列归一化
%% 得到观测向量y（已优化for循环）
                y=zeros(sampling_num,1); %储存观测变量y
                y_i=1:sampling_num;
                y(y_i)=a(sample(y_i)+slip-sampling_range); %得到观测向量y
%% 调用SAMP方法得到方程系数
                theta=SAMP(y,A,1); %调用SAMP方法得到方程系数
%% 逆归一化，得到方程系数（已优化for循环）
                coefficient=zeros(1,coefficient_num); %储存方程系数
                target_samp_i=1:coefficient_num;
                coefficient(target_samp_i)=theta(target_samp_i)./l(target_samp_i); %逆归一化，得到方程系数
%% 储存每次循环的所有系数矩阵及采样点
                all_target(i,:)=coefficient; %储存每次循环的所有系数矩阵
                all_sample(i,:)=sample; %储存每次循环的所有采样点
%% 得到重构数据（已优化for循环）
                a_rebuild_resolution=zeros(sampling_range,index_sum); %储存代数并乘系数后的测量矩阵
                h=1:sampling_range; %循环得到滑动窗口内的重构数据
%% 得到代数并乘系数后的测量矩阵
                for a_rebuild_i=1:index_sum
                    re2=1;
                    for jjj =2:varn  
                        re2=re2.*(var{jjj}(h+slip-sampling_range).^(index_num(a_rebuild_i,jjj-1)));
                    end
                    a_rebuild_resolution(:,a_rebuild_i)=coefficient(1,a_rebuild_i)*re2;
%                     a_rebuild_resolution(:,a_rebuild_i)=coefficient(1,a_rebuild_i)*(b(h+slip-sampling_range).^(index_num(a_rebuild_i,1))).*(c(h+slip-sampling_range).^(index_num(a_rebuild_i,2))).*(d(h+slip-sampling_range).^(index_num(a_rebuild_i,3)));
                end
                a_rebuild=sum(a_rebuild_resolution,2); %储存重构数据
%% 得到平均重构相对误差
                a_real=a((slip-sampling_range+1):(slip),1); %得到采样范围内各个点的实际数据
                a_error=(a_rebuild-a_real)./a_real; %得到采样范围内各个点的重构相对误差（由于窗口内的数据确认是无异常的，所以分母为实际数据）
                all_mean_error(i,1)=sum(abs(a_error))/sampling_range; %储存所有平均重构相对误差，平均重构相对误差等于所有重构相对误差绝对值之和的平均值
                if sum(abs(a_error)<mean_error_limit)==(sampling_range)
                    all_a_error(i,1)=1;
                end
            end
%% 找到所有平均重构相对误差符合要求的方程位置及个数
            all_mean_error_pass=find(all_a_error==1); %找到所有平均重构相对误差符合要求的方程位置
            all_mean_error_pass_num=length(all_mean_error_pass); %找到所有平均重构相对误差符合要求的方程个数
        end
%% 选择平均重构相对误差最小的若干个方程
        [result,index]=sort(all_mean_error); %平均重构相对误差排序
        all_mean_error_pass_part=index(1:choose_num,1); %选择方程采集数量（根据平均重构相对误差由小到大选择）
%% 得到所有经过选择的方程在测试点的修复数据以及修复相对误差
        test_num=slip; %设置测试点位置
        test_error=zeros(choose_num,1); %储存测试点修复相对误差
        for r=1:choose_num %循环求解所有经过选择的方程在测试点的修复数据
        target=all_target(all_mean_error_pass_part(r,1),:); %得到方程系数
        test_rebuild_01=solve_equationsYUAN(var_bcd,target,test_num,index_sum,index_num); %得到测试点的修复数据
        test_real_01=a(test_num,1); %得到测试数据的实际数据
        test_error_01=abs((test_rebuild_01-test_real_01)/test_real_01); %储存测试点的修复相对误差（由于测试点的数据确认是无异常的，所以分母为实际数据）
        test_error(r,1)=test_error_01; %储存测试点的修复相对误差（即测试数据的修复数据与实际数据的相对误差的绝对值）
        end
%% 储存所有的测试点修复相对误差
        all_test_error(:,slip-slip_start+1)=test_error; %储存所有的测试点修复相对误差
%% 在满足while条件下得到下一时刻的修复数据
        a_correct=0; %设置a_correct初始值，进入循环
        test_error_judge=0; %设置test_error_judge初始值，进入循环
        while(a_correct==0&&test_error_judge<test_error_judge_limit) %若下一时刻的修复数据不为零或由小到大排列的测试点修复相对误差大于等于阈值，则跳出循环
            minposition=find(test_error==min(test_error)); %找到当前最小的测试点修复相对误差的位置
            only_minposition=minposition(1,1); %得到当前唯一的最小修复相对误差的位置
            test_error_judge=test_error(only_minposition); %得到当前最小的测试点修复相对误差
            test_error(minposition)=1; %当前的最小测试点修复相对误差置1
            if test_error_judge<test_error_judge_limit  %重要判断，若当前最小修复相对误差小于阈值，修复下一时刻的数据
                target=all_target(all_mean_error_pass_part(only_minposition,1),:); %得到方程系数
                a_correct=solve_equationsYUAN(var_bcd,target,test_num+1,index_sum,index_num); %得到下一时刻的修复数据
            end
        end
%% 储存所有的当前最小修复相对误差的位置
        all_only_minposition(slip-slip_start+1,1)=only_minposition; %储存所有的当前最小修复相对误差的位置
%% 防止循环次数溢出
        loop_num=loop_num-1; %下一次进入判断平均重构相对误差符合要求的方程数量是否足够的循环时，保证循环次数不变（上一次循环已经证明当前循环次数足够）
%% 判断修复数据是否符合要求，若不符合要求则将所有修复数据储存起来，以备筛选
%% 重要判断，若下一时刻的修复数据不为0，则进入以下流程
        if a_correct~=0
            relative_error_judge=abs((a_correct-a(test_num+1))/a_correct); %下一时刻的修复相对误差（由于下一点的数据不确定是否异常的，所以分母为修复数据）
%% 若相对误差大于阈值，且判断次数小于阈值，则进入以下流程（在判断次数的阈值内得到小于阈值的修复相对误差，说明修复数据满足要求）
            if relative_error_judge>=relative_error_judge_limit&&judge_times<judge_num 
%% 储存备选的修复数据以及修复相对误差
                alternative_a_correct(judge_times+1,1)=a_correct; %储存备选的修复数据
                alternative_relative_error_judge(judge_times+1,1)=relative_error_judge; %储存备选的修复相对误差
%% 储存备选的系数矩阵及采样点
                alternative_final_target(judge_times+1,:)=target; %储存备选的系数矩阵
%% 跳出循环，重构进入循环
                judge_times=judge_times+1; %判断次数+1
                relative_error_judge_limit=relative_error_judge_limit*error_expansion_coefficient; %缩小下一点修复相对误差限制（*修复对误差扩张系数）
                test_error_judge=1; %为continue重新进入循环设置初值
                continue; %重新进入判断测试点修复相对误差是否符合要求的循环中
            end
        end
    end
%% 数据修复过程（1）
    if judge_times<judge_num %若判断次数小于阈值，说明在规定次数内得到了符合要求的修复数据
        a(slip+1)=a_correct; %直接将修复数据代入，完成修复
% a(slip+1)=a(test_num+1); %数据无异常，不用修复
%% 储存最终的系数矩阵
        final_target(1+slip-slip_start,:)=target; %储存最终的系数矩阵
%% 储存第一次出现的方程及位置（在final_target中）
        first_equations_combination(1+slip-slip_start,:)=target; %储存第一次出现的方程（在final_target中）
        first_equations_combination_local(1+slip-slip_start,1)=1+slip-slip_start; %储存第一次出现的方程位置（在final_target中）
    elseif judge_times==judge_num  %若判断次数等于阈值，说明在规定次数内没有得到符合要求的修复数据，进入筛选流程
        [j_result,j_index]=sort(alternative_relative_error_judge); %对修复数据相对误差进行排序
%% 选取中间点（无论异常是否确定，都将修复数据代入，完成修复）
        if sum(alternative_relative_error_judge>=alternative_relative_error_judge_limit)==judge_num
            middle_num=ceil(judge_num/2); %若所有可选修复相对误差大于或等于阈值的数量等于总数量，说明监测数据很可能出现异常，中间点取中间值
        elseif sum(alternative_relative_error_judge>=alternative_relative_error_judge_limit)<judge_num
            middle_num=1; %若所有可选修复相对误差大于或等于阈值的数量小于总数量，说明监测数据可能无异常，造成相对误差较大的原因可能是数据噪声，中间点取最小值    
        end
%% 数据修复过程（2）
        find_min_num=j_index(middle_num,1); %选择排第几小的修复相对误差，得到其位置
        a_correct=alternative_a_correct(find_min_num,1); %得到所选的修复数据
        a(slip+1)=a_correct; %将修复数据代入，完成修复
%% 储存最终的系数矩阵
        final_target(1+slip-slip_start,:)=alternative_final_target(find_min_num,:); %储存最终的系数矩阵
%% 储存第一次出现的方程及位置（在final_target中）
        first_equations_combination(1+slip-slip_start,:)=alternative_final_target(find_min_num,:); %储存第一次出现的方程
        first_equations_combination_local(1+slip-slip_start,1)=1+slip-slip_start; %储存第一次出现的方程位置（在final_target）
    end
%% 储存当前备选的修复数据
    all_alternative_a_correct(:,1+slip-slip_start)=alternative_a_correct; % 储存当前备选的修复数据
%% 储存当前备选的修复相对误差
    all_alternative_relative_error_judge(:,1+slip-slip_start)=alternative_relative_error_judge; % 储存当前备选的修复相对误差
%% 储存当前循环次数及判断次数
    loop_num=loop_num+1; %恢复循环次数
    all_loop_num(1+slip-slip_start,1)=loop_num; %储存当前循环次数
    all_judge_times(1+slip-slip_start,1)=judge_times; %储存当前判断次数
end
%% 进度条关闭
close(hwait); %进度条关闭%% 实际数据和修复数据处理
%% 实际数据和修复数据处理
x_final(1,:)=aa((slip_start+1):(slip_end+1)); %得到实际数据 
x_final(2,:)=a((slip_start+1):(slip_end+1)); %得到修复数据
normal_repair(1,:)=x_final(2,:)-x_final(1,:); %得到实际数据与修复数据的误差
relative_normal_repair(1,:)=normal_repair(1,:)./x_final(1,:); %得到实际数据与修复数据的相对误差（用于检验算法，已知实际数据都是无异常的，看修复数据精度如何）
relative_normal_repair(2,:)=normal_repair(1,:)./x_final(2,:); %得到实际数据与修复数据的相对误差（由于修复数据更可靠，所以分母为修复数据）
%% 实际数据与修复数据的对比图
figure; %画图开始
plot(1:rebuild_num,x_final(1,:),'ks-','MarkerSize',3,'MarkerFaceColor','k'); %画出实际数据曲线
hold on; %图像合并
plot(1:rebuild_num,x_final(2,:),'r^-','MarkerSize',3,'MarkerFaceColor','r'); %画出修复数据曲线
% axis([1 rebuild_num 1.10 1.40]); %设置坐标轴范围
legend('实际数据','修复数据'); %图像标注
hold off; %画图结束
%% 修复相对误差曲线
figure; %画图开始
plot(1:rebuild_num,relative_normal_repair(1,:),'ko-','MarkerSize',3,'MarkerFaceColor','k'); %画出修复相对误差曲线
hold off; %画图结束
%% 结束计时
toc; %结束计时