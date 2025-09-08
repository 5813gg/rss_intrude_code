% 清空工作区变量、命令窗口，并关闭所有图形窗口
clear;    % 清除工作区所有变量
clc;      % 清空命令窗口
close all; % 关闭所有图形窗口

% 设置基本参数
Fn = 200;        % 采样频率200Hz
set_mw = 0.1;    % RSSI处理的最小值阈值
T = 8;           % 采样时间8秒
num_samples = Fn * T;  % 计算总采样点数
time_axis = linspace(0, T, num_samples); % 生成时间轴

% 加载静态环境RSSI数据
[~, ~, ~, rssis] = load_data("D:\code\data\1201\static.dat");
 
% 读取CSI（信道状态信息）数据
csi_trace = read_bf_file("D:\code\data\1201\static.dat");

% 初始化AGC（自动增益控制）数组
agc = zeros(1, size(csi_trace, 2));

% 提取每帧的AGC值
for ii = 1:size(agc, 2)
    agc(:, ii) = csi_trace{ii}.agc;
end

% RSSI数据处理
rssi_mag = 0;  % 初始化RSSI幅度
rssi_mag = rssi_mag + dbinv(rssis(:, 1).'); % 将dB转换为线性单位
rss = db(rssi_mag, 'pow') - 44 - agc; % 计算真实的RSS值（减去固定偏移和AGC影响）

% 提取中间稳定的RSS数据段（去除前后不稳定部分）
RSS0 = rss(Fn+1:9*Fn);

% 初始化去噪后的静态数据变量
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

% 归一化滤波器系数
Rayleightemp = Rayleightemp / sum(Rayleightemp);
% 找到滤波器系数的最大值及其位置
[maxr, max_position] = max(Rayleightemp);

 
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

% FFT分析参数设置
Fs = 200;       % 采样频率200Hz
windowSize = 400; % 窗口大小400点
RSS0length = length(rssi0_smooth); % 静态数据长度
N = windowSize;  % FFT点数

% 初始化能量比数组
energy_below_30Hz = zeros(1,RSS0length-windowSize);

% 滑动窗口计算10Hz以下能量占比
for i = 1:RSS0length-windowSize
    % 执行FFT变换
    Y = fft(rssi0_smooth(i:i+windowSize-1));
    % 计算频率轴
    f = (0:N-1)*(Fs/N);
    % 计算双边谱幅度
    P2 = abs(Y/N);
    % 转换为单边谱
    P1 = P2(1:N/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    % 生成频率轴
    f1 = f(1:N/2+1);
    f1 = f1(2:end); % 移除DC分量
    P1 = P1(2:end); % 相应地移除DC分量幅度
    
    % 计算能量
    energy_total = sum(P1); % 总能量
    indices = f1 < 10;      % 10Hz以下频率索引
    energy_below_10Hz = sum(P1(indices)); % 10Hz以下能量
    
    % 计算10Hz以下能量占比
    percentage_below_10Hz(1,i) = (energy_below_10Hz / energy_total);
end

% 获取最大能量占比作为检测阈值
[threshold, max_idx] = max(percentage_below_10Hz);

% 加载不同距离的动态RSS数据（走廊环境tx-rx1）
RSS1 = get_rss("D:\code\rssi\rss_data\corridor\txrx1\1.dat");
RSS2 = get_rss("D:\code\rssi\rss_data\corridor\txrx1\2.dat");
RSS3 = get_rss("D:\code\rssi\rss_data\corridor\txrx1\3.dat");
RSS4 = get_rss("D:\code\rssi\rss_data\corridor\txrx1\4.dat");
RSS5 = get_rss("D:\code\rssi\rss_data\corridor\txrx1\5.dat");
RSS6 = get_rss("D:\code\rssi\rss_data\corridor\txrx1\6.dat");
RSS7 = get_rss("D:\code\rssi\rss_data\corridor\txrx1\7.dat");
RSS8 = get_rss("D:\code\rssi\rss_data\corridor\txrx1\8.dat");
RSS9 = get_rss("D:\code\rssi\rss_data\corridor\txrx1\9.dat");
RSS10 = get_rss("D:\code\rssi\rss_data\corridor\txrx1\10.dat");



% 组织动态RSS数据（1-8米）
RSS_all = {RSS1,RSS2,RSS3,RSS4,RSS5,RSS6,RSS7,RSS8};

% 初始化去噪后的动态数据变量
RSS_denoised_mw = cell(size(RSS_all));
RSS_denoised_all = cell(size(RSS_all));
RSS_denoised_all_mw = cell(size(RSS_all));

% 动态RSS数据处理（与静态数据处理类似）
for ii = 1:size(RSS_denoised_mw, 2)
    for jj = 1:size(RSS_all{ii}, 2)
        RSS_denoised_mw{ii}(:, jj) = 10^(RSS_all{ii}(:, jj)/10);
    end
end

for ii = 1:size(RSS_denoised_all, 2)
    for jj = 1:size(RSS_denoised_mw{ii}, 2)
        RSS_denoised_all_mw{ii}(:, jj) = (RSS_denoised_mw{ii}(:, jj)-min(RSS_denoised_mw{ii}(RSS_denoised_mw{ii}~=0)))/min(RSS_denoised_mw{ii}(RSS_denoised_mw{ii}~=0));
    end
end

for ii = 1:size(RSS_denoised_all, 2)
    RSS_denoised_all_mw{ii}(RSS_denoised_all_mw{ii}==0) = min(RSS_denoised_all_mw{ii}(RSS_denoised_all_mw{ii}~=0))-set_mw;
end

for ii = 1:size(RSS_denoised_all, 2)
    for jj = 1:size(RSS_denoised_all_mw{ii}, 2)
        RSS_denoised_all{ii}(:, jj) = 10*log10(RSS_denoised_all_mw{ii}(:, jj));
    end
end

% 使用瑞利滤波器平滑动态RSS数据
rssi_smooth = cell(size(RSS_denoised_all));
for ii = 1:size(rssi_smooth, 2)
    for jj = 1:size(RSS_denoised_all{ii}, 2)
        if jj < max_position
            rssi_smooth{ii}(:, jj) = [zeros(1, max_position-jj), RSS_denoised_all{ii}(:, 1:jj+2*r-1-max_position)]*Rayleightemp';
        elseif jj+2*r-1-max_position > size(RSS_denoised_all{ii}, 2)
            rssi_smooth{ii}(:, jj) = [RSS_denoised_all{ii}(:, jj-max_position+1:size(RSS_denoised_all{ii}, 2)), zeros(1, jj+2*r-1-max_position-size(RSS_denoised_all{ii}, 2))]*Rayleightemp';
        else
            rssi_smooth{ii}(:, jj) = RSS_denoised_all{ii}(:, jj-max_position+1 : jj+2*r-1-max_position)*Rayleightemp';
        end
    end
end

% 初始化检测结果数组
RSS_detection = zeros(1, size(rssi_smooth, 2));

% 对每个距离的动态数据进行检测
for ii = 1:size(RSS_detection, 2)
    RSS_intrusion = rssi_smooth{ii};
    RSS_intrusion_length = length(RSS_intrusion);
    
    % 滑动窗口计算10Hz以下能量占比
    for i = 1:RSS_intrusion_length-windowSize
        Y = fft(rssi_smooth{ii}(i:i+windowSize-1));
        f = (0:N-1)*(Fs/N);
        P2 = abs(Y/N);
        P1 = P2(1:N/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f1 = f(1:N/2+1);
        f1 = f1(2:end);
        P1 = P1(2:end);
        
        % 计算能量比
        energy_total = sum(P1);
        indices = f1 < 10;
        energy_below_10Hz = sum(P1(indices));
        percentage_below_10Hz_intrusion(1,i) = (energy_below_10Hz / energy_total);
    end
    
    % 计算检测率（超过阈值的比例）
    num_greater_than_threshold = sum(percentage_below_10Hz_intrusion > threshold);
    detection_percent = num_greater_than_threshold/length(percentage_below_10Hz_intrusion(~isnan(percentage_below_10Hz_intrusion)));
    RSS_detection(1,ii) = detection_percent*100;
end

% 绘制检测结果图
figure;
hold on;
plot(RSS_detection, 'b-o', 'LineWidth', 2);
xticklabels({'1','2','3','4','5','6','7','8','9','10'});  
set(gcf, 'Color', 'white');
set(gca, 'Color', 'white');
xlabel('距离(米)');
ylabel('检测率(%)');
box(gca, 'on');