function [csi,amplitude,phase,rssi] = load_data(path)
% Input ： Path=存储dat文件的路径
% Output： 输出csi，amplitude，phase和rssi的值
%          csi相关数据维度为 [接收天线*发送天线*30个子载波数据维度  *  时间长度] ==》270*T
% Example：csi,amplitude,phase,rssi = load("Data\11-23\user1_01.dat")
    csi_trace=read_bf_file(path);
    L=length(csi_trace);

    csi = zeros(L,270); 
    rssi = zeros(L,3);
    perm = zeros(L,3);
    
    for m=1:L
        csi_entry=csi_trace{m};
        csi_data=get_scaled_csi(csi_entry); 
        %csi_data = getfield (csi_entry, 'csi'); 
        
        % CSI数据录入
        % tx=1,rx=1:3
        csi(m,1:30)=csi_data(1,1,:);
        csi(m,31:60)=csi_data(1,2,:);
        csi(m,61:90)=csi_data(1,3,:);
        % tx=2,rx=1:3
        csi(m,91:120)=csi_data(2,1,:);
        csi(m,121:150)=csi_data(2,2,:);
        csi(m,151:180)=csi_data(2,3,:);
        % tx=3,rx=1:3
        csi(m,181:210)=csi_data(3,1,:);
        csi(m,211:240)=csi_data(3,2,:);
        csi(m,241:270)=csi_data(3,3,:);
        
        % RSSI数据录入
        rssi(m,1) = csi_entry.rssi_a;
        rssi(m,2) = csi_entry.rssi_b;
        rssi(m,3) = csi_entry.rssi_c;
    end
    
    amplitude = abs(csi);
    phase = angle(csi);
end
% csi_data=[];
% csi_label = [];
% 
% csi_data[i]= amp..[i]'=;
% label[i]= 
% 
% end
% 
% save(csi_fall,csi_daata,label);



