function [ output_args ] = scatter_gauss( file )
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
g = load(file);
scatter(g(1,4), g(1, 5));
hold on
scatter(g(2:27, 4), g(2:27, 5));
hold on 
scatter(g(28:63, 4), g(28:63, 5));
hold on
scatter(g(64:112, 4), g(64:112, 5));
hold on
scatter(g(113:176, 4), g(113:176, 5));
hold on

end

