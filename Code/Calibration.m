clear;clc;close all;
N = 18;     % 摆放位置的个数，共有N个摆放位置
M = 88;     % 角点数量
img_width=680; img_height=700; 
id_normal=[1,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,21];
%%%%%%%%%%%%%%%%%%%%数据读取
%1)加载相机内参 和 每幅相移图的RT矩阵
S = load ('..\Data\cameraParams.mat');
K=S.cameraParams.IntrinsicMatrix;K=K';  %读取相机内参数
fx = K(1,1);  fy = K(2,2);  u0 = K(1,3);  v0 = K(2,3); s_factor=K(1,2);
k1=S.cameraParams.RadialDistortion;	% 相机畸变系数
k2=S.cameraParams.TangentialDistortion;
k=[k1,k2];
R=zeros(3,3,N);T=zeros(1,3,N); 
for i=1:N
    R(:,:,i)=S.cameraParams.RotationMatrices(:,:,i); 
    T(:,:,i)=S.cameraParams.TranslationVectors(i,:);  
end
%2)加载角点世界坐标
World_coordinate = [S.cameraParams.WorldPoints, zeros(M,1)];
%3)加载相位图
img_phase = zeros(N,img_width,img_height); 
for i=1:N
    temp=load(['..\Data\Phase\',num2str(id_normal(i)),'.mat']);
    img_phase(i,:,:) = temp.phi_unwrapped';     
end
%%%%%%%%%%%%%%%%%%%%求解N个平面的方程
Camera_coordinate=zeros(M,3,N);
plane=zeros(N,4);
for i=1:N
    for j=1:M
        Camera_coordinate(j,:,i)=World_coordinate(j,:)*R(:,:,i)+T(:,:,i);  %求取相机坐标
    end
    D = [ Camera_coordinate(:,1,i), Camera_coordinate(:,2,i), Camera_coordinate(:,3,i),ones(size(Camera_coordinate(:,1,i)))];
    [U,S,V] = svd(D,0);
    plane(i,:) = V(:,end)';  % 平面系数 [A B C D]
end
%%%%%%%%%%%%%%%%%%%%求解直线方程
line=zeros(2,img_width,img_height);
for v=1:img_height
    for u=1:img_width
        ydn=(v-v0)/fy;
        xdn=(u-u0-s_factor*ydn)/fx;
        r=sqrt(xdn^2+ydn^2);
        xn=xdn-xdn*(k(1)*r^2+k(2)*r^4);   
        yn=ydn-ydn*(k(1)*r^2+k(2)*r^4);
        zn=1;
        line(:,u,v)=[xn/zn;yn/zn];   
    end
end
%%%%%%%%%%%%%%%%%%%%直线与平面的交点，共有N*img_width*img_height个交点(Xc,Yc,Zc)
Camera_coordinate=zeros(N,3,img_width,img_height);
for v=1:img_height
    for u=1:img_width
        Camera_coordinate(:,3,u,v)=-plane(:,4)./(plane(:,1)*line(1,u,v)+plane(:,2)*line(2,u,v)+plane(:,3));
        Camera_coordinate(:,1,u,v)=line(1,u,v)*Camera_coordinate(:,3,u,v);
        Camera_coordinate(:,2,u,v)=line(2,u,v)*Camera_coordinate(:,3,u,v); 
    end
end
%%%%%%%%%%%%%%%%%%%%逐像素标定
Calibration_Table=zeros(6,img_width,img_height);
for v=1:img_height
    for u=1:img_width
        phases_uv=img_phase(:,u,v);    %uv对应的N个相位
        Camera_coordinates_uv=Camera_coordinate(:,:,u,v);    %N*3
        %标定Zc中的参数a1,a2,a3,a4
        D = [ phases_uv, ones(size(phases_uv)), -phases_uv.*Camera_coordinates_uv(:,3), -Camera_coordinates_uv(:,3)];
        [U,S,V] = svd(D,0); 
        Calibration_Table(1:4,u,v) = V(:,end); 
        %标定Xc中的参数a5
        D = [ Camera_coordinates_uv(:,1), Camera_coordinates_uv(:,3)];
        [U,S,V] = svd(D,0);temp=V(:,end); 
        Calibration_Table(5,u,v)=-temp(2)/temp(1);
        %标定Yc中的参数a6
        D = [Camera_coordinates_uv(:,2), Camera_coordinates_uv(:,3)];
        [U,S,V] = svd(D,0);temp=V(:,end); 
        Calibration_Table(6,u,v)=-temp(2)/temp(1);
    end
end
%%%%%%%%%%%%%%%%%%%%标定参数表保存
save('..\Data\Calibration_Table.mat', 'Calibration_Table');

%%%%%%%%%%%%%%%%%%%%像素误差分析
S=load ('..\Data\Calibration_Table.mat');
Calibration_Table=S.Calibration_Table;
pixel_wise_error = zeros(img_width, img_height);
for v=1:img_height
    for u=1:img_width
        calculate_position_z = (img_phase(:,u,v)*Calibration_Table(1,u,v)+Calibration_Table(2,u,v))./(img_phase(:,u,v)*Calibration_Table(3,u,v)+Calibration_Table(4,u,v));
        calculate_position_x = Calibration_Table(5,u,v)*calculate_position_z;
        calculate_position_y = Calibration_Table(6,u,v)*calculate_position_z;      
        x_dev = calculate_position_x-Camera_coordinate(:,1,u,v);
        y_dev = calculate_position_y-Camera_coordinate(:,2,u,v);
        z_dev = calculate_position_z-Camera_coordinate(:,3,u,v);
        pixel_wise_error(u,v) = sum(sqrt(z_dev.^2))/N;
    end
end
error_vector = reshape(pixel_wise_error,[1,img_height*img_width]);
error_vector(error_vector>0.3)=0.3;
mean_error = mean(error_vector);
max_error=max(error_vector);
std_error = std(error_vector);

%%%%%%%%%%%%%%%%%%%%画误差的频率分布直方图
bargram = histogram(error_vector,100);bargram.FaceColor='b';bargram.EdgeColor=[0.5 0.5 0.5];%频数直方图
xlim([0 0.17]);
%%%%%%%%%%%%%%%%%%%%像素分级
levelNum=5; %最多分类5个ABCDE
Label_Table = zeros(img_width, img_height);
frequency=bargram.Values; errors=bargram.BinEdges;
mode_error=errors(find(frequency==max(frequency))+1);
Label_Table(pixel_wise_error<mode_error-std_error)=1;  %第一类
Label_Table(pixel_wise_error>=mode_error+(levelNum-3)*std_error)=levelNum;  %最后一类
for i=2:levelNum-1
  Label_Table((pixel_wise_error>=mode_error+(i-3)*std_error)&(pixel_wise_error<mode_error+(i-2)*std_error))=i;  %最后一类
end
save('..\Data\Label_Table.mat', 'Label_Table');