f=2.4e9;%载波频率
c=3e8;%光速
lambda=c/f;%波长
Pt=0.02;%发射功率为20dB // ESP32 AP发射功率范围（10-20dBm）可通过编程调整
n=2;%路径损耗因子
Gt=3;Gr=5;%发射增益，接收增益
l=0.4 ; w=0.2 ; h=1.7 ;%人体长宽高
ht=1;%发射机与地面的距离
Tbl=0.8;Tmo=0.8;Tf=0.8;%设置反射损耗系数
% x=1:0.1:3;
% y=2*x-1;
%x=2:0.1:3;
y=1:0.1:5;
x = ones(1,length(y)) * 2;


 x0=1;y0=1;%人的初始位置
% x=1:0.1:3.5;
%y=1:0.1:2.5;
TX=[0,0];
RX=[4,0];
device={TX,RX};
Ht=0;%初始化
Hs=0;Hbl=0;Hmo=0;
Gl=Gt*Gr;
RSSt=zeros(1,size(x,2));
RSS0t=zeros(1,size(x,2));
sblt=zeros(1,size(x,2));
smot=zeros(1,size(x,2));



b=6;%发射机到南墙的距离
a=7;%发射机到东墙的距离
kw=b/a;%发射器到右下墙的斜率
angle=90;
for tt=1:length(x)
    vertex1=[x(tt)+(w/2)*sin(angle)-(l/2)*cos(angle),y(tt)-(l/2)*sin(angle)-(w/2)*cos(angle)];
    vertex2=[x(tt)+(w/2)*sin(angle)+(l/2)*cos(angle),y(tt)+(l/2)*sin(angle)-(w/2)*cos(angle)];
    vertex3=[x(tt)-(w/2)*sin(angle)+(l/2)*cos(angle),y(tt)+(l/2)*sin(angle)+(w/2)*cos(angle)];
    vertex4=[x(tt)-(w/2)*sin(angle)-(l/2)*cos(angle),y(tt)-(l/2)*sin(angle)+(w/2)*cos(angle)];
    rectangle={vertex1,vertex2,vertex3,vertex4};
    area=[l*h,w*h,l*h,w*h];
    %人体四面中点的坐标
    center1=(vertex1+vertex2)/2;
    center2=(vertex2+vertex3)/2;
    center3=(vertex3+vertex4)/2;
    center4=(vertex4+vertex1)/2;
    center={center1,center2,center3,center4};
    %四点的斜率
    k=[vertex1(2)/vertex1(1),vertex2(2)/vertex2(1),vertex3(2)/vertex3(1),vertex4(2)/vertex4(1)];
    kmin=min(k);
    kmax=max(k);
    Sbl1=max(1/max(kmin,kw)-1/kmax,0)*(b^2*(h-ht)/y0+ht);%投影在南墙的面积
    Sbl2=max(min(kmax,kw)-kmin,0)*(a^2*(h-ht)/x0+ht);%投影在东墙的面积
    Sbl=Sbl1+Sbl2;%阻挡面积
    sblt(:,tt)=Sbl;
    
    score=zeros(1,size(center,2));
    dot=0;
    scorecenter=0;
    scoretotal=0;
    %判断反射面
    reflectsurface={vertex2-vertex1,vertex3-vertex2,vertex4-vertex3,vertex1-vertex4};
    for ii=1:size(center,2)
        dot=0;
        for jj=1:size(reflectsurface,2)
            % 跳过当前行
            if jj == ii
                continue;
            end
            for kk=1:size(device,2)
                txcenter=center{ii}-device{kk};
                rectangletx=device{kk}-rectangle{jj};
                rectanglecenter=center{ii}-rectangle{jj};
                txrectangle=rectangle{jj}-device{kk};
                txrectangle_other=reflectsurface{jj}-rectangletx;
                A=reflectsurface{jj}(1)*rectangletx(2)-reflectsurface{jj}(2)*rectangletx(1);
                B=reflectsurface{jj}(1)*rectanglecenter(2)-reflectsurface{jj}(2)*rectanglecenter(1);
                result1=A*B;
                %
                C=txcenter(1)*txrectangle(2)-txcenter(2)*txrectangle(1);
                D=txcenter(1)*txrectangle_other(2)-txcenter(2)*txrectangle_other(1);
                result2=C*D;
                if result1<=0 && result2<=0
                    dot=1;
                else
                    dot=0;
                end
                scorecenter=dot+scorecenter;
            end
            scoretotal=scorecenter+scoretotal;

        end
        score(:,ii)=scoretotal;
        scoretotal=0;
        scorecenter=0;
    end
    position=find(score == 0);
    Smo=area(position);
    smot(:,tt)=Smo;

    Ht=0;%初始化
    Hs=0;Hbl=0;Hmo=0;
    Gl=Gt*Gr;
    dbl=sqrt((x(tt)-TX(1))^2+(y(tt)-TX(2))^2);
    dmo=sqrt((x(tt)-TX(1))^2+(y(tt)-TX(2))^2)+sqrt((x(tt)-RX(1))^2+(y(tt)-RX(2))^2);

    abl=lambda./(4*pi)*(1/dbl^(n/2));
    amo=lambda./(4*pi)*(1/dmo^(n/2));
    Hbl=abl*Sbl*Tbl*exp(-1i*2*pi*dbl/lambda);
    Hmo=amo*Smo*Tmo*exp(-1i*2*pi*dmo/lambda);
    Hs=-0.012185782739616 + 0.024635206843988i;
    zzzz=abs(Hs)^2*Pt;
    Ht=Hs-Hbl+Hmo;
    RSS0=10*log10(abs(Hs)^2*Pt);
    RSS=10*log10(abs(Ht)^2*Pt);
    RSSt(:,tt)=RSS;
    RSS0t(:,tt)=RSS0;
end







