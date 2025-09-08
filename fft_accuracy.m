% 清空工作区变量、命令窗口，并关闭所有图形窗口
clear;    % 清除工作区变量
clc;      % 清空命令窗口
close all; % 关闭所有图形窗口

% 设置基本参数
Fn = 200;        % 采样频率（Hz）
set_mw = 0.1;    % RSSI处理的最小值阈值
T = 1/Fn;        % 采样周期（秒）

% 从数据文件加载静态RSSI数据
[~, ~, ~, rssis] = load_data("D:\code\data\1201\static.dat");

% 读取CSI（信道状态信息）数据文件
csi_trace = read_bf_file("D:\code\data\1201\static.dat");

% 初始化AGC（自动增益控制）数组
agc = zeros(1, size(csi_trace, 2));

% 提取每帧的AGC值
for ii = 1:size(agc, 2)
    agc(:, ii) = csi_trace{ii}.agc;
end

% 初始化RSSI幅度并转换单位
rssi_mag = 0;  % 初始化RSSI幅度
rssi_mag = rssi_mag + dbinv(rssis(:, 1).'); % 将dB转换为线性单位

% 计算真实的RSS值（减去固定偏移和AGC影响）
rss = db(rssi_mag, 'pow') - 44 - agc;

% 提取中间稳定的RSS数据段（去除前后不稳定部分）
RSS0 = rss(Fn+1:9*Fn);

% 初始化结果存储矩阵
zzz = zeros(5, 10); 

% 加载不同距离的动态RSS数据文件（1-10米）
RSS1_ALL = files('D:\code\data\FFT_DATA\1\*.dat');
RSS2_ALL = files('D:\code\data\FFT_DATA\2\*.dat');
RSS3_ALL = files('D:\code\data\FFT_DATA\3\*.dat');
RSS4_ALL = files('D:\code\data\FFT_DATA\4\*.dat');
RSS5_ALL = files('D:\code\data\FFT_DATA\5\*.dat');
RSS6_ALL = files('D:\code\data\FFT_DATA\6\*.dat');
RSS7_ALL = files('D:\code\data\FFT_DATA\7\*.dat');
RSS8_ALL = files('D:\code\data\FFT_DATA\8\*.dat');
RSS9_ALL = files('D:\code\data\FFT_DATA\9\*.dat');
RSS10_ALL = files('D:\code\data\FFT_DATA\10\*.dat');

% 以下是注释掉的志愿者数据路径（备用）
% RSS1_ALL=files('D:\code\data\1227\volunteer\zsy\1\*.dat');
% ...

% 初始化静态数据处理变量
RSS0_denoised = zeros(size(RSS0));  % 去噪后的静态数据
RSS0_mw = zeros(size(RSS0));        % 转换为毫瓦单位的静态数据
RSS0_denoised_mw = zeros(size(RSS0_mw)); % 归一化后的静态数据

% 将静态RSS从dB转换为线性单位（毫瓦）
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
    % 瑞利分布公式计算滤波器系数
    Rayleightemp(i) = (i-1)/(r_sigma^2) * exp(-(i-1)^2/(2*r_sigma^2));
end

% 绘制滤波器系数
plot(Rayleightemp)
% 归一化滤波器系数
Rayleightemp = Rayleightemp / sum(Rayleightemp);
% 找到滤波器系数的最大值及其位置
[maxr, max_position] = max(Rayleightemp);

% 使用瑞利滤波器平滑静态RSS数据
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

% 绘制平滑后的静态数据
plot(rssi0_smooth)

% 设置FFT分析参数
Fs = 200; % 采样频率为200Hz
windowSize = length(rssi0_smooth); % 窗口大小
RSS0length = length(rssi0_smooth); % 静态数据总长度
N = length(rssi0_smooth);

% 对静态数据进行FFT分析
Y = fft(rssi0_smooth); % 执行FFT变换
f = (0:N-1)*(Fs/N);    % 计算频率轴
P2 = abs(Y/N);         % 归一化幅度
P1 = P2(1:N/2+1);      % 取单边谱
P1(2:end-1) = 2*P1(2:end-1); % 处理双边谱
f1 = f(1:N/2+1);       % 相应的频率范围
f1 = f1(2:end);        % 移除DC分量（0 Hz）
P1 = P1(2:end);        % 相应地移除DC分量幅度

% 计算频谱能量
energy_total = sum(P1); % 总能量

% 计算10Hz以下频段的能量
indices = f1 < 10; 
energy_below_10Hz = sum(P1(indices)); 

% 计算10Hz以下能量占比（作为后续检测阈值）
threshold = (energy_below_10Hz / energy_total);

% 将所有距离的动态数据放入cell数组
RSSall = {RSS1_ALL, RSS2_ALL, RSS3_ALL, RSS4_ALL, RSS5_ALL, RSS6_ALL, RSS7_ALL, RSS8_ALL, RSS9_ALL, RSS10_ALL};

% 初始化精度数组
accuracy_all_f = zeros(1, size(RSSall, 2));
accuracy_all_f_n = zeros(1, size(RSSall, 2));

% 处理每个距离的动态数据
for zz = 1:size(RSSall, 2)
    RSS_all = RSSall{zz};
    
    % 初始化去噪后的数据存储cell
    RSS_denoised_all = cell(size(RSS_all));
    RSS_denoised_mw = cell(size(RSS_all));
    RSS_denoised_all_mw = cell(size(RSS_denoised_mw));
    
    % 将动态RSS从dB转换为线性单位（毫瓦）
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
    
    % 使用瑞利滤波器平滑每个距离的动态RSS数据
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
    
    % 初始化能量比数组
    sum_ratio = zeros(1, size(rssi_smooth, 2));
    
    % 对每个动态数据样本进行FFT分析
    for ii = 1:size(rssi_smooth, 2)
        RSS_intrusion = rssi_smooth{ii};
        RSS_intrusion_length = length(RSS_intrusion); 
        
        % 执行FFT变换
        Y = fft(RSS_intrusion); 
        f = (0:RSS_intrusion_length-1)*(Fs/RSS_intrusion_length); 
        P2 = abs(Y/RSS_intrusion_length);
        P1 = P2(1:RSS_intrusion_length/2+1); 
        P1(2:end-1) = 2*P1(2:end-1); 
        f1 = f(1:RSS_intrusion_length/2+1); 
        f1 = f1(2:end); 
        P1 = P1(2:end); 
        
        % 计算能量比
        energy_total = sum(P1);
        indices = f1 < 10; 
        energy_below_10Hz = sum(P1(indices)); 
        percentage_below_10Hz_intrusion = (energy_below_10Hz / energy_total);
        
        % 存储结果
        sum_ratio(:, ii) = percentage_below_10Hz_intrusion;
    end
    
    % 将当前距离的结果存入矩阵
    zzz(:, zz) = sum_ratio';
end

% 准备绘图数据
test = zzz;

% 绘制箱线图（1-9米的数据）
boxplot(test(:, 1:end-1));

% 绘制静态数据的参考线
therods = ones(1, 9)*threshold;
hold on
plot(therods, 'r--')

% 设置图形属性
xlabel('距离(米)');
ylabel('10Hz以下能量占比');
set(findobj(gca, 'type', 'line'), 'linewidth', 1.5) % 设置线宽
set(gcf, 'Color', 'white');  % 设置图形背景为白色
set(gca, 'Color', 'white');  % 设置坐标轴背景为白色