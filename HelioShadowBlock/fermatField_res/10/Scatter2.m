function [ output_args ] = Scatter2( file_name )
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
data = load(file_name);
hpx = data(:,1);
hpy = data(:,2);
sd_bk = data(:,3);
flux = data(:,4);
rela_dis = data(:,5);
approx_dis = data(:,6);

subplot(2,2,1);
scatter(hpx, hpy,20, flux, 'filled');
subplot(2,2,2);
scatter(hpx, hpy, 20, sd_bk, 'filled');
subplot(2,2,3);
scatter(hpx, hpy, 20, rela_dis, 'filled');
subplot(2,2,4);
scatter(hpx, hpy, 20, approx_dis, 'filled');

end

