function [ output_args ] = flux( ~ )
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
grid = load('grid.txt');
gauss = load('gauss_51.txt');
grid_recv = load('grid_51.txt');



x1 = grid(:,1);
y1=grid(:,2);
z1=grid(:,3);
x2=gauss(:,1);
y2=gauss(:,2);
z2=gauss(:,3);
x3=grid_recv(:,1);
y3=grid_recv(:,2);
z3=grid_recv(:,3);

subplot(1,2,1);

scatter(x2,y2,10, z2,'filled');
caxis([-0.2,1.7]);
%subplot(1,3,2);
%scatter(x2,y2,20,z2,'filled');
subplot(1,2,2);

scatter(x3,y3,10,z3,'filled');
caxis([-0.2,1.7]);

end

