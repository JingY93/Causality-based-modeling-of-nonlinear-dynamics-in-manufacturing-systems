%% 函数定义
% 设有m个盒子,每个盒子里各有n1,n2,...,nm个不同的球。
% 从每个盒子中各取一个球，共有几种组合，列出所有组合。
function COM = Combination_MN(num_ball)
%% 参数设置
% num_ball 每个盒子中球的个数，[n1,n2,...,nm];
num_box = length(num_ball); % 盒子个数
%% 函数体
NumCom = 1; % 共有多少种组合
for i_box = 1:num_box
    NumCom = NumCom*num_ball(i_box);
end
COM = zeros(NumCom,num_box); % 初始化组合向量
i_com = 1; % 当前组合在COM中的地址
com = ones(1,num_box); % 当前组合
COM(i_com,:) = com; % 所有组合
while ~isequal(com,num_ball)
    com(end) = com(end) + 1;
    for i_box = num_box:-1:2
        if com(i_box)>num_ball(i_box)
            com(i_box) = 1; com(i_box-1) = com(i_box-1) + 1;
        end
    end
    i_com = i_com + 1;
    COM(i_com,:) = com;
end
COM=COM-1;
end

