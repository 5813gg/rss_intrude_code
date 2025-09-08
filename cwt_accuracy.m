% 清空工作区变量、命令窗口，并关闭所有图形窗口
clear;    % 清除工作区变量
clc;      % 清空命令窗口
close all; % 关闭所有图形窗口

% 设置基本参数
Fn = 200;        % 采样频率（Hz）
windowsize = 200; % 处理窗口大小
set_mw = 0.1;    % RSSI处理的最小值阈值

% 从数据文件加载RSSI数据
[~, ~, ~, rssis] = load_data("D:\code\rssi\rss_data\cwt\CWT_DATA\static.dat");

% 读取CSI（信道状态信息）数据文件
csi_trace = read_bf_file("D:\code\rssi\rss_data\cwt\CWT_DATA\static.dat");

% 初始化AGC（自动增益控制）数组
agc = zeros(1, size(csi_trace, 2));

% 提取每帧的AGC值
for ii = 1:size(agc, 2)
    agc(:, ii) = csi_trace{ii}.agc;
end

% 初始化RSSI幅度并转换单位
rssi_mag = 0;
rssi_mag = rssi_mag + dbinv(rssis(:, 1).'); % 将dB转换为线性单位

% 计算真实的RSS值（减去固定偏移和AGC影响）
rss = db(rssi_mag, 'pow') - 44 - agc;

% 提取中间稳定的RSS数据段（去除前后不稳定部分）
RSS0 = rss(Fn+1:9*Fn);

% 加载不同距离的RSS数据文件（1-10米）
RSS1_ALL = files('D:\code\rssi\rss_data\cwt\CWT_DATA\1\*.dat');
RSS2_ALL = files('D:\code\rssi\rss_data\cwt\CWT_DATA\2\*.dat');
RSS3_ALL = files('D:\code\rssi\rss_data\cwt\CWT_DATA\3\*.dat');
RSS4_ALL = files('D:\code\rssi\rss_data\cwt\CWT_DATA\4\*.dat');
RSS5_ALL = files('D:\code\rssi\rss_data\cwt\CWT_DATA\5\*.dat');
RSS6_ALL = files('D:\code\rssi\rss_data\cwt\CWT_DATA\6\*.dat');
RSS7_ALL = files('D:\code\rssi\rss_data\cwt\CWT_DATA\7\*.dat');
RSS8_ALL = files('D:\code\rssi\rss_data\cwt\CWT_DATA\8\*.dat');
RSS9_ALL = files('D:\code\rssi\rss_data\cwt\CWT_DATA\9\*.dat');
RSS10_ALL = files('D:\code\rssi\rss_data\cwt\CWT_DATA\10\*.dat');

% 初始化结果存储矩阵
zzz = zeros(5, 10);

% 初始化去噪后的RSS数组
RSS0_denoised = zeros(size(RSS0));
RSS0_mw = zeros(size(RSS0));
RSS0_denoised_mw = zeros(size(RSS0_mw));

% 将RSS从dB转换为线性单位（毫瓦）
for ii = 1:size(RSS0_mw, 2)
    RSS0_mw(:, ii) = 10^(RSS0(:, ii)/10);
end

% 归一化处理（减去最小值后除以最小值）
for ii = 1:size(RSS0_denoised_mw, 2)
    RSS0_denoised_mw(:, ii) = (RSS0_mw(:, ii)-min(RSS0_mw(RSS0_mw~=0)))/min(RSS0_mw(RSS0_mw~=0));
end

% 处理零值，用最小值减去阈值替代
RSS0_denoised_mw(RSS0_denoised_mw==0) = min(RSS0_denoised_mw(RSS0_denoised_mw~=0))-set_mw;

% 将处理后的数据转换回dB单位
for ii = 1:size(RSS0_denoised, 2)
    RSS0_denoised(:, ii) = 10*log10(RSS0_denoised_mw(:, ii));
end

% 设置瑞利滤波参数
r = 9;          % 滤波器半宽
r_sigma = 4;    % 瑞利分布参数

% 生成瑞利滤波器系数
Rayleightemp = ones(1, r*2-1);
for i = 1:r*2-1
    % 瑞利分布公式
    Rayleightemp(i) = (i-1)/(r_sigma^2) * exp(-(i-1)^2/(2*r_sigma^2));
end

% 绘制滤波器系数
plot(Rayleightemp)
% 归一化滤波器系数
Rayleightemp = Rayleightemp / sum(Rayleightemp);
% 找到滤波器系数的最大值及其位置
[maxr, max_position] = max(Rayleightemp);

% 使用瑞利滤波器平滑RSS数据
rssi0_smooth = zeros(size(RSS0_denoised));
for ii = 1:size(rssi0_smooth, 2)
    % 处理边界情况
    if ii < max_position
        % 左边界：前面补零
        rssi0_smooth(:, ii) = [zeros(1, max_position-ii), RSS0_denoised(:, 1:ii+2*r-1-max_position)]*Rayleightemp';
    elseif ii+2*r-1-max_position > size(rssi0_smooth, 2)
        % 右边界：后面补零
        rssi0_smooth(:, ii) = [RSS0_denoised(:, ii-max_position+1:size(rssi0_smooth, 2)), zeros(1, ii+2*r-1-max_position-size(rssi0_smooth, 2))]*Rayleightemp';
    else
        % 正常情况：直接应用滤波器
        rssi0_smooth(:, ii) = RSS0_denoised(:, ii-max_position+1 : ii+2*r-1-max_position)*Rayleightemp';
    end
end

% 设置小波分解参数
wname = 'db4';  % 使用db4小波
level = 5;      % 分解层数

% 对平滑后的RSS进行小波分解
[rssi0_C, rssi0_L] = wavedec(rssi0_smooth, level, wname);

% 初始化小波系数阈值处理
c_thresh_rss0 = {0, 0, 0};
% 提取近似系数
rssi0_A3 = appcoef(rssi0_C, rssi0_L, wname, level);

% 对各层细节系数进行阈值处理
for i = 1:level
    % 获取当前层的细节系数
    cD_rss0 = detcoef(rssi0_C, rssi0_L, i);
    % 使用rigrsure方法计算阈值
    thr_rss0 = thselect(cD_rss0, 'rigrsure');
    % 应用软阈值处理
    cD_rss0 = wthresh(cD_rss0, 's', thr_rss0);
    % 存储处理后的系数
    c_thresh_rss0{i} = cD_rss0;
end

% 重组小波系数（近似系数+各层处理后的细节系数）
cl_rss0 = [rssi0_A3 c_thresh_rss0{level} c_thresh_rss0{level-1} c_thresh_rss0{level-2} c_thresh_rss0{level-3} c_thresh_rss0{level-4}];
% 小波重构
rssi0_dwt = waverec(cl_rss0, rssi0_L, wname);

% 设置连续小波变换参数
wavename = 'db4';   % 使用db4小波
totalscal = 100;    % 总尺度数
% 计算小波的中心频率
wcf = centfrq(wavename);
% 计算尺度参数
cparam = 2*wcf*totalscal;
scal = cparam./(1:totalscal);
% 将尺度转换为频率
freq = scal2frq(scal, wavename, 1/Fn);

% 对原始RSS进行连续小波变换
coefs_rss0_n = cwt(RSS0, scal, wavename);
% 计算每个尺度的最大系数
sp_rss0_n = max(abs(coefs_rss0_n), [], 2);
% 找到最大系数对应的尺度和位置
[sp_max_rss0_n, sp_position_rss0_n] = max(sp_rss0_n);
% 获取对应的频率
freq_interest_rss0_n = freq(sp_position_rss0_n);
% 计算小波分解层数
N_rss0_n = ceil(log2((Fn/2)/freq_interest_rss0_n));
% 进行小波分解
[rss0_c_n, rss0_l_n] = wavedec(RSS0, N_rss0_n, wname);
% 提取细节系数
detail_coefficient_rss0_n = detcoef(rss0_c_n, rss0_l_n, N_rss0_n);
% 计算细节系数的方差
var_detail0_n = var(detail_coefficient_rss0_n);

% 对去噪后的RSS进行连续小波变换
coefs_rss0 = cwt(rssi0_dwt, scal, wavename);
% 计算每个尺度的最大系数
sp_rss0 = max(abs(coefs_rss0), [], 2);
% 找到最大系数对应的尺度和位置
[sp_max_rss0, sp_position_rss0] = max(sp_rss0);
% 获取对应的频率
freq_interest_rss0 = freq(sp_position_rss0);
% 计算小波分解层数
N_rss0 = ceil(log2((Fn/2)/freq_interest_rss0));
% 进行小波分解
[rss0_c, rss0_l] = wavedec(rssi0_dwt, N_rss0, wname);
% 提取细节系数
detail_coefficient_rss0 = detcoef(rss0_c, rss0_l, N_rss0);
% 计算细节系数的方差
var_detail0 = var(detail_coefficient_rss0);

% 将所有距离的数据放入cell数组
RSSall = {RSS1_ALL, RSS2_ALL, RSS3_ALL, RSS4_ALL, RSS5_ALL, RSS6_ALL, RSS7_ALL, RSS8_ALL, RSS9_ALL, RSS10_ALL};

% 初始化精度数组
accuracy_all_c = zeros(1, size(RSSall, 2));
accuracy_all_c_n = zeros(1, size(RSSall, 2));

% 处理每个距离的数据
for zz = 1:size(RSSall, 2)
    RSS_all = RSSall{zz};
    
    % 初始化去噪后的数据存储cell
    RSS_denoised_all = cell(size(RSS_all));
    RSS_denoised_mw = cell(size(RSS_all));
    RSS_denoised_all_mw = cell(size(RSS_denoised_mw));
    
    % 将RSS从dB转换为线性单位（毫瓦）
    for ii = 1:size(RSS_denoised_mw, 2)
        for jj = 1:size(RSS_all{ii}, 2)
            RSS_denoised_mw{ii}(:, jj) = 10^(RSS_all{ii}(:, jj)/10);
        end
    end
    
    % 归一化处理（减去最小值后除以最小值）
    for ii = 1:size(RSS_denoised_all, 2)
        for jj = 1:size(RSS_denoised_mw{ii}, 2)
            RSS_denoised_all_mw{ii}(:, jj) = (RSS_denoised_mw{ii}(:, jj)-min(RSS_denoised_mw{ii}(RSS_denoised_mw{ii}~=0)))/min(RSS_denoised_mw{ii}(RSS_denoised_mw{ii}~=0));
        end
    end
    
    % 处理零值，用最小值减去阈值替代
    for ii = 1:size(RSS_denoised_all, 2)
        RSS_denoised_all_mw{ii}(RSS_denoised_all_mw{ii}==0) = min(RSS_denoised_all_mw{ii}(RSS_denoised_all_mw{ii}~=0))-set_mw;
    end
    
    % 将处理后的数据转换回dB单位
    for ii = 1:size(RSS_denoised_all, 2)
        for jj = 1:size(RSS_denoised_all_mw{ii}, 2)
            RSS_denoised_all{ii}(:, jj) = 10*log10(RSS_denoised_all_mw{ii}(:, jj));
        end
    end
    
    % 使用瑞利滤波器平滑每个距离的RSS数据
    rssi_smooth = cell(size(RSS_denoised_all));
    for ii = 1:size(rssi_smooth, 2)
        for jj = 1:size(RSS_denoised_all{ii}, 2)
            % 处理边界情况
            if jj < max_position
                % 左边界：前面补零
                rssi_smooth{ii}(:, jj) = [zeros(1, max_position-jj), RSS_denoised_all{ii}(:, 1:jj+2*r-1-max_position)]*Rayleightemp';
            elseif jj+2*r-1-max_position > size(RSS_denoised_all{ii}, 2)
                % 右边界：后面补零
                rssi_smooth{ii}(:, jj) = [RSS_denoised_all{ii}(:, jj-max_position+1:size(RSS_denoised_all{ii}, 2)), zeros(1, jj+2*r-1-max_position-size(RSS_denoised_all{ii}, 2))]*Rayleightemp';
            else
                % 正常情况：直接应用滤波器
                rssi_smooth{ii}(:, jj) = RSS_denoised_all{ii}(:, jj-max_position+1 : jj+2*r-1-max_position)*Rayleightemp';
            end
        end
    end
    
    % 对平滑后的数据进行小波去噪
    rssi_dwt = cell(size(rssi_smooth));
    for ii = 1:size(rssi_dwt, 2)
        % 小波分解
        [C, L] = wavedec(rssi_smooth{ii}, level, wname);
        
        % 初始化小波系数阈值处理
        c_thresh_rss = {0, 0, 0};
        % 提取近似系数
        A3 = appcoef(C, L, wname, level);
        
        % 对各层细节系数进行阈值处理
        for jj = 1:level
            % 获取当前层的细节系数
            cD_rss = detcoef(C, L, jj);
            % 使用rigrsure方法计算阈值
            thr_rss = thselect(cD_rss, 'rigrsure');
            % 应用软阈值处理
            cD_rss = wthresh(cD_rss, 's', thr_rss);
            % 存储处理后的系数
            c_thresh_rss{jj} = cD_rss;
        end
        
        % 重组小波系数（近似系数+各层处理后的细节系数）
        cl_rss = [A3 c_thresh_rss{level} c_thresh_rss{level-1} c_thresh_rss{level-2} c_thresh_rss{level-3} c_thresh_rss{level-4}];
        % 小波重构
        rssi_dwt{ii} = waverec(cl_rss, L, wname);
    end
    
    % 计算原始RSS数据的细节系数方差
    var_detail_n = zeros(1, size(RSS_all, 2));
    for ii = 1:size(var_detail_n, 2)
        % 连续小波变换
        coefs_rss = cwt(RSS_all{ii}, scal, wavename);
        % 计算每个尺度的最大系数
        sp_rss = max(abs(coefs_rss), [], 2);
        % 找到最大系数对应的尺度和位置
        [sp_max_rss, sp_position_rss] = max(sp_rss);
        % 获取对应的频率
        freq_interest_rss = freq(sp_position_rss);
        % 计算小波分解层数
        N_rss = ceil(log2((Fn/2)/freq_interest_rss));
        % 进行小波分解
        [rss_c, rss_l] = wavedec(RSS_all{ii}, N_rss, wname);
        % 提取细节系数
        detail_coefficient_rss = detcoef(rss_c, rss_l, N_rss);
        % 计算细节系数的方差
        var_detail_n(:, ii) = var(detail_coefficient_rss);
    end
    
    % 计算去噪后RSS数据的细节系数方差
    var_detail = zeros(1, size(rssi_dwt, 2));
    for ii = 1:size(var_detail, 2)
        % 连续小波变换
        coefs_rss = cwt(rssi_dwt{ii}, scal, wavename);
        % 计算每个尺度的最大系数
        sp_rss = max(abs(coefs_rss), [], 2);
        % 找到最大系数对应的尺度和位置
        [sp_max_rss, sp_position_rss] = max(sp_rss);
        % 获取对应的频率
        freq_interest_rss = freq(sp_position_rss);
        % 计算小波分解层数
        N_rss = ceil(log2((Fn/2)/freq_interest_rss));
        % 进行小波分解
        [rss_c, rss_l] = wavedec(rssi_dwt{ii}, N_rss, wname);
        % 提取细节系数
        detail_coefficient_rss = detcoef(rss_c, rss_l, N_rss);
        % 计算细节系数的方差
        var_detail(:, ii) = var(detail_coefficient_rss);
    end
    
    % 存储结果（这里选择存储原始数据的方差）
    zzz(:, zz) = var_detail_n';
end

% 绘制箱线图（1-9米的数据）
boxplot(zzz(:, 1:end-1));
hold on
% 绘制参考线（使用原始RSS的细节系数方差）
therods = ones(1, 9)*var_detail0_n';
plot(therods, 'r--')

% 设置图形属性
xlabel('距离(米)');
ylabel('细节系数方差');
xticklabels({'1', '2', '3', '4', '5', '6', '7', '8', '9'}); 
set(gcf, 'Color', 'white');  % 设置图形背景为白色
set(gca, 'Color', 'white');  % 设置坐标轴背景为白色