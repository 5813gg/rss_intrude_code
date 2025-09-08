% 清空工作区变量、命令窗口，并关闭所有图形窗口
clear;    % 清除工作区变量
clc;      % 清空命令窗口
close all; % 关闭所有图形窗口

% 设置基本参数
Fn = 200;        % 采样频率（Hz）
set_mw = 0.1;    % RSSI处理的最小值阈值

% 从数据文件加载静态RSSI数据
[~, ~, ~, rssis] = load_data("D:\code\rssi\rss_data\fft_max\static.dat");

% 读取CSI（信道状态信息）数据文件
csi_trace = read_bf_file("D:\code\rssi\rss_data\fft_max\static.dat");

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

% 加载不同距离的动态RSS数据（1-10米）
RSS1 = get_rss("D:\code\rssi\rss_data\fft_max\1.dat");
RSS2 = get_rss("D:\code\rssi\rss_data\fft_max\2.dat");
RSS3 = get_rss("D:\code\rssi\rss_data\fft_max\3.dat");
RSS4 = get_rss("D:\code\rssi\rss_data\fft_max\4.dat");
RSS5 = get_rss("D:\code\rssi\rss_data\fft_max\5.dat");
RSS6 = get_rss("D:\code\rssi\rss_data\fft_max\6.dat");
RSS7 = get_rss("D:\code\rssi\rss_data\fft_max\7.dat");
RSS8 = get_rss("D:\code\rssi\rss_data\fft_max\8.dat");
RSS9 = get_rss("D:\code\rssi\rss_data\fft_max\9.dat");
RSS10 = get_rss("D:\code\rssi\rss_data\fft_max\10.dat");

% 定义距离边界（10米到1米）
boud = [10,9,8,7,6,5,4,3,2,1];

% 将RSS数据按距离从远到近存储（10m-1m）
RSS_all = {RSS10,RSS9,RSS8,RSS7,RSS6,RSS5,RSS4,RSS3,RSS2,RSS1};

% 计算采样周期
T = 1/Fn;

% 初始化FFT幅度数组
RSS_amplitude = zeros(1, size(RSS_all, 2));
RSSamplitude = cell(size(RSS_all));

% 对每个距离的RSS数据进行FFT分析
for ii = 1:size(RSS_amplitude, 2)
    L = length(RSS_all{ii});  % 获取数据长度
    t = (0:L-1)*T;           % 生成时间轴
    
    % 执行FFT变换
    RSS_fft = fft(RSS_all{ii});
    
    % 计算双边谱幅度
    P2 = abs(RSS_fft/L);
    
    % 转换为单边谱
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    % 生成频率轴
    f = Fn*(0:(L/2))/L;
    
    % 去除直流分量（0Hz）
    P1(1) = 0;
    
    % 存储FFT结果
    RSSamplitude{ii} = P1;
    
    % 记录最大幅度值
    RSS_amplitude(:, ii) = max(P1);
end

% 对静态RSS数据进行FFT分析
L0 = length(RSS0);        % 静态数据长度
t0 = (0:L0-1)*T;          % 时间轴
RSS0_fft = fft(RSS0);     % FFT变换

% 计算静态数据的FFT幅度
P20 = abs(RSS0_fft/L0);   % 双边谱
P10 = P20(1:L0/2+1);      % 单边谱
P10(2:end-1) = 2*P10(2:end-1); % 幅度修正

% 生成频率轴
f0 = Fn*(0:(L0/2))/L0;

% 去除直流分量
P10(1) = 0;

% 获取静态数据的最大FFT幅度及其位置
[RSS0_amplitude, index0] = max(P10);

% 找出第一个超过静态阈值的距离位置
position = find(RSS_amplitude > RSS0_amplitude, 1, 'first');

% 计算可检测的最远距离
dmax = boud(:, position);

% 输出检测结果
fprintf("能够检测到人员入侵的最远距离为：%.2f米", dmax);

% 生成参考线数据（静态数据的FFT幅度）
therods = ones(size(RSS_amplitude(:, 1:end))) * RSS0_amplitude;

% 绘制结果图形
figure;
hold on;

% 绘制各距离的FFT最大幅度（按1-10米顺序）
plot(flip(RSS_amplitude), 'b-o', 'LineWidth', 2);

% 绘制静态数据参考线
plot(therods, 'r--', 'LineWidth', 2);

% 设置X轴标签（距离）
xticklabels({'1','2','3','4','5','6','7','8','9','10'});  

% 设置图形背景为白色
set(gcf, 'Color', 'white');
set(gca, 'Color', 'white');

% 添加坐标轴标签
xlabel('Distance(m)');
ylabel('FFT amplitude');

% 显示坐标轴边框
box(gca, 'on');