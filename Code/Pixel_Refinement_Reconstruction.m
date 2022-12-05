%% 三维重建
clc; close all; clear;
img_width=680; img_height=700; 
N=9;
S=load ('..\Data\Calibration_Table.mat');
Calibration_Table=S.Calibration_Table;
Label_Table_Struct = load('..\Data\Label_Table.mat');
Label_Table = Label_Table_Struct.Label_Table;
%%%%%%%%%%%%%%%%%%%%%%相位图读取
Label='ABCDE';
for k=1:5    %标签的等级总数
    for idx = 9   %Selct the position_id
        temp=load(['..\Data\Phase\Sphere\',num2str(idx),'.mat']);
        img_phase= temp.phi_unwrapped';     %注意左侧索引顺序为(u,v),此处有转置
%         imshow(img_phase',[]);
        %%%%%%%%%%%%%%%%%%%%重建 
        Mask_Table_Struct = load(['..\Data\Phase\Sphere\Mask_',num2str(idx),'.mat']);
        Mask=Mask_Table_Struct.Mask_Phase;
        data=zeros(img_height*img_width,3);
        for v=1:img_height
            for u=1:img_width
                %%%%%通过判断label来选择是否重建
                if Label_Table(u,v)<=k & Mask(u,v)==1
                    calculate_position_z = (img_phase(u,v)*Calibration_Table(1,u,v)+Calibration_Table(2,u,v))./(img_phase(u,v)*Calibration_Table(3,u,v)+Calibration_Table(4,u,v));
                    calculate_position_x = Calibration_Table(5,u,v)*calculate_position_z;
                    calculate_position_y = Calibration_Table(6,u,v)*calculate_position_z; 
                    data(img_width*(v-1)+u,:)=[calculate_position_x,calculate_position_y,calculate_position_z];
                else
                    data(img_width*(v-1)+u,:)=nan;
                end
            end
        end
        data1=data(all(~isnan(data),2),:); %去掉为nan的数据
        ptCloud = pointCloud(data1);
        figure; pcshow(ptCloud);
    end
end