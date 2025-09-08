 
% 函数功能：从CSI数据文件中提取并处理RSS（接收信号强度）数据
% 输入参数：file_name - 数据文件名（.dat文件）
% 输出参数：RSS - 处理后的RSS信号（1600个采样点）
function RSS=get_rss(file_name)
Fn=200; % 设置采样频率为200Hz（用于后续数据截取）

% 调用load_data函数加载数据（忽略前三个输出，只获取rssi数据）

[~,~,~,rssis]=load_data(file_name);% 表示忽略该输出参数，rssis获取三个接收天线的RSSI值
csi_trace=read_bf_file(file_name); % 直接读取CSI原始数据文件
agc=zeros(1,size(csi_trace,1)); % 初始化AGC（自动增益控制）数组
% 大小为1×数据包数量（每个数据包一个AGC值）

% 遍历所有数据包，获取每个包的AGC值
for ii=1:size(agc,2)
    agc(:,ii)=csi_trace{ii}.agc; % 从每个CSI数据包中提取AGC值
end
 rssi_mag=0; % 初始化RSSI幅度变量（线性单位）
 % 将RSSI从dBm转换为线性幅度（使用第3个接收天线的数据）
 % 注释掉的代码是使用第1个天线的版本
 %rssi_mag = rssi_mag + dbinv(rssis(:,1).');
 rssi_mag = rssi_mag + dbinv(rssis(:,3).'); % 使用天线3的RSSI
 % 计算最终的RSS值（转换为dBm单位）：
% 1. db(...,'pow')将线性功率转换回dBm
% 2. -44是Intel 5300网卡的固定偏移量
% 3. -agc减去自动增益控制值
 rss = db(rssi_mag, 'pow') - 44 - agc;

 % 截取中间稳定的1600个采样点（去掉开头200和结尾200的不稳定数据）
% Fn=200, 所以 Fn+1:9*Fn 就是 201:1800
 RSS=rss(Fn+1:9*Fn);
end

