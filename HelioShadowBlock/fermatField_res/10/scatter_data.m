function [ output_args ] = scatter_data( file_name )
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
file = load(file_name);
x = file(:,1);
y = file(:,2);
sd = file(:,3);
sum = file(:,4);

scatter(x, y, 20, sum, 'filled');

end

