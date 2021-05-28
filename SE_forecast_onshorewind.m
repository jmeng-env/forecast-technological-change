% copyright Xiangping and Jing, 2020

clear;clc;close all;
year=xlsread('D:\Nexus365\Rupert Way - forecast_comparison_paper\R2 revision\data and figures\Fig 3 data.xlsx','Onshore wind', 'A2:A38');
year=floor(year); % ignore the month
price=xlsread('D:\Nexus365\Rupert Way - forecast_comparison_paper\R2 revision\data and figures\Fig 3 data.xlsx','Onshore wind', 'C2:C38');
logprice=log(price);
sl=diff(logprice);
z=xlsread('D:\Nexus365\Rupert Way - forecast_comparison_paper\R2 revision\data and figures\Fig 3 data.xlsx','Onshore wind', 'B2:B38');
logz = log(z);
startyear=1983;
endyear=2019;
m=endyear-startyear+1;

%% fit the Normal
% mu and sigma
pd_sl= fitdist(sl, 'normal');
histfit(sl, 10, 'normal')
%% predict

year_pred = 2020:2030;
delta_t = year_pred - year(end);

% theoretical mu and std
mu_pred = pd_sl.mu * delta_t;
sd_pred = sqrt(pd_sl.sigma^2 * delta_t);

% sample
n = 10000;
sl_mc = normrnd(pd_sl.mu, pd_sl.sigma, [n,length(year_pred)]);
sl_pred = cumsum(sl_mc,2);
price_pred = price(end) * exp(sl_pred);


%% figure: percentiles

price_pred_median = median(price_pred);
price_pred_05th = prctile(price_pred,05);
price_pred_10th = prctile(price_pred,10);
price_pred_50th = prctile(price_pred,50);
price_pred_90th = prctile(price_pred,90);
price_pred_95th = prctile(price_pred,95);
year_all = [year;year_pred'];
price_all_median = [price;price_pred_median'];


figure
semilogy(year_all,price_all_median,'-')
hold on
semilogy(year, price,'bo')
semilogy(year_pred, price_pred_median,'rx')
semilogy(year_pred, price_pred_05th,'r')
semilogy(year_pred, price_pred_10th,'r')
semilogy(year_pred, price_pred_50th,'r')
semilogy(year_pred, price_pred_90th,'r')
semilogy(year_pred, price_pred_95th,'r')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% wright
%% fit the Normal to change of log(z)
% logcapacity=logz;
% 
% year2=[year,ones(m,1)];
% [bbbb,bint2,r2,rint2,regress2]=regress(logcapacity(23:32),year2(23:32,:));


%%

% X01=0:(2030-endyear);
% X01=1:(2019-endyear);
% X02=year2(:,1);
% 
% x=[X02,ones(m,1)];
% ycapacity=[logcapacity];
% 
% % prediction to 2050
% yCalc1_his = x*bbbb;
% yCalc1= bbbb(1,:)*X01+ ycapacity(end);

% sample and then predict z
 ratio=(z(37)/z(32))^(1/5)-1;
z_pred_median = zeros(2030-2019,1);
for i=1:size(z_pred_median)
    z_pred_median(i)=z(end)*(1+ratio)^(i);
end


%% fit the normao to change of logprice/logz

n=10000;
sl2 = diff(logprice)./diff(logz);
pd_sl2 = fitdist(sl2, 'normal');

% sample and then predict price
sl2_mc = normrnd(pd_sl2.mu, pd_sl2.sigma, [n,length(z_pred_median)]);
price_pred_w = zeros(n, length(z_pred_median));
for i = 1:length(z_pred_median)
    if i == 1
        price_pred_w(:,i) = price(end) .* (z_pred_median(i)/z(end)) .^ sl2_mc(:,i);
    else
        price_pred_w(:,i) = price_pred_w(:,i-1) .* (z_pred_median(i)./z_pred_median(i-1)) .^ sl2_mc(:,i);
    end
end

%% figure: percentiles
% price_pred=exp(price_pred);
price_pred_median_w = median(price_pred_w);
price_pred_05th_w = prctile(price_pred_w,05);
price_pred_10th_w = prctile(price_pred_w,10);
price_pred_50th_w = prctile(price_pred_w,50);
price_pred_90th_w = prctile(price_pred_w,90);
price_pred_95th_w = prctile(price_pred_w,95);
year_all = [year;year_pred'];
price_all_median = [price;price_pred_median_w'];

% figure
% plot(year_all,price_all_median,'-')
% hold on
% plot(year, price,'bo')
% plot(year_pred, price_pred_median,'rx')
% plot(year_pred, price_pred_05th,'r-.')
% plot(year_pred, price_pred_25th,'r-.')
% plot(year_pred, price_pred_75th,'r-.')
% plot(year_pred, price_pred_95th,'r-.')
% hold off

figure

Y_W=xlsread('D:\Nexus365\Rupert Way - forecast_comparison_paper\R2 revision\Fig3_W1_M1_2030_forecast_data\Fig_3_PDW_Onshore_wind.xlsx', 'C2:I12'); % global $/w
y_W=zeros(7,12);
y_W(:,1)=exp(logprice(end));
y_W(:,2:12)=Y_W';

Y_M=xlsread('D:\Nexus365\Rupert Way - forecast_comparison_paper\R2 revision\Fig3_W1_M1_2030_forecast_data\Fig_3_PDM_Onshore_wind.xlsx', 'C2:I12'); % global $/w
y_M=zeros(7,12);
y_M(:,1)=exp(logprice(end));
y_M(:,2:12)=Y_M';
ee=[1.258585608	1.521128469	1.723925886	1.944011931	2.291998016];

semilogy(year_all,price_all_median,'-')
hold on
semilogy(year, price,'bo')
semilogy(year_pred, price_pred_median_w,'bx')
semilogy(year_pred, price_pred_05th_w,'b')
semilogy([year(end) year_pred(1)], [price(end) price_pred_05th_w(1)],'b')
semilogy(year_pred, price_pred_10th_w,'b')
semilogy([year(end) year_pred(1)], [price(end) price_pred_10th_w(1)],'b')
semilogy(year_pred, price_pred_50th_w,'b')
semilogy([year(end) year_pred(1)], [price(end) price_pred_50th_w(1)],'b')
semilogy(year_pred, price_pred_90th_w,'b')
semilogy([year(end) year_pred(1)], [price(end) price_pred_90th_w(1)],'b')
semilogy(year_pred, price_pred_95th_w,'b')
semilogy([year(end) year_pred(1)], [price(end) price_pred_95th_w(1)],'b')

semilogy(year_pred, price_pred_median,'gx')
semilogy(year_pred, price_pred_05th,'g')
semilogy([year(end) year_pred(1)], [price(end) price_pred_05th(1)],'g')
semilogy(year_pred, price_pred_10th,'g')
semilogy([year(end) year_pred(1)], [price(end) price_pred_10th(1)],'g')
semilogy(year_pred, price_pred_50th,'g')
semilogy([year(end) year_pred(1)], [price(end) price_pred_50th(1)],'g')
semilogy(year_pred, price_pred_90th,'g')
semilogy([year(end) year_pred(1)], [price(end) price_pred_90th(1)],'g')
semilogy(year_pred, price_pred_95th,'g')
semilogy([year(end) year_pred(1)], [price(end) price_pred_95th(1)],'g')
semilogy(endyear:2030, y_W, 'r-') %%W1
semilogy(endyear:2030, y_M, 'k-') %%W1
semilogy(2030, ee, 'p-') %%W1
data=zeros(5,4);
data(:,1)=y_W(2:6,12);
data(:,2)=y_M(2:6,12);
data(:,3)=prctile(price_pred_w(:,11),[5 25 50 75 95]);
data(:,4)=prctile(price_pred(:,11),[5 25 50 75 95]);

hold on 

ylim([0.1, 10])
xlim([1970, 2030])
hold off

 
