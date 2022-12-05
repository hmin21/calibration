clc; close all; clear;
img_width=680; img_height=700; 
%%%%%%%%%%%%%%%%%%%%%%相位图读取
id_normal=22:26;  %复杂曲面图像
S=load ('..\Data\Calibration_Table.mat');
Calibration_Table=S.Calibration_Table;
for id=1:5   % Select the reconstruction target
    temp=load(['..\Data\Phase\',num2str(id_normal(id)),'.mat']);
    img_phase= temp.phi_unwrapped';     %注意左侧索引顺序为(u,v),此处有转置
    figure;imshow(img_phase',[]);colorbar;caxis([20,70]);
    %%%%%%%%%%%%%%%%%%%%重建
    data=zeros(img_height*img_width,3);
    for v=1:img_height
        for u=1:img_width
            calculate_position_z = (img_phase(u,v)*Calibration_Table(1,u,v)+Calibration_Table(2,u,v))./(img_phase(u,v)*Calibration_Table(3,u,v)+Calibration_Table(4,u,v));
            calculate_position_x = Calibration_Table(5,u,v)*calculate_position_z;
            calculate_position_y = Calibration_Table(6,u,v)*calculate_position_z; 
            data(img_width*(v-1)+u,:)=[calculate_position_x,calculate_position_y,calculate_position_z];
        end
    end
    data (data(:,3) <165,:) = nan;
    data (data(:,3) >185,:) = nan;
    ptCloud = pointCloud(data);
    % pcwrite(data, '../Reconstruction/球/ptcloud.pcd', 'Encoding', 'ascii'); %将程序中的xyz数据写入pcd文件中
    % pc = pcread('../Reconstruction/球/ptcloud.pcd');
    figure; pcshow(ptCloud); axis off; 
end