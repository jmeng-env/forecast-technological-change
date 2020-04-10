 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% copyright Jing, 2018
 %%%  1, change data source
 %%%  2, change n, startyear, endyear
 %%%  3, check the distribution
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;close all;
filepath='D:\OneDrive - University Of Cambridge\Cambridge\INNOPATH\learning curve data\summary\';
wind=xlsread([filepath,'Wind.xlsx'], 'installed cost', 'B7:D41'); % global $/w
year=wind(:,1);
logyear=year;
price=wind(:,3);
logprice=log(price);

startyear=1983;
endyear=2017;
n=endyear-startyear+1;

year2=[year,ones(n,1)];
[bbbb,bint2,r2,rint2,regress2]=regress(logprice,year2);

% total years lag3
t=5;
predictor=cell(n-t,n-t);
predict=cell(n-t,n-t);
for j=1:n-t
    for i=1:n-(t-1)-j
      predict{j,i}=logprice(i:n+1-j,:);                         %e
      predictor{j,i}=[logyear(i:n+1-j),ones((n+2)-j-i,1)];  
    end
end

result={};
b1=zeros((n-(t-1))*(n-t)/2,1);
m=0;
for i=1:n-t
  for j=1:n-(t-1)-i
    [b,bint,r,rint,stats]=regress(predict{i,j},predictor{i,j});
    result{i,j,1}=b;
    result{i,j,2}=bint;
    result{i,j,3}=r;
    result{i,j,4}=rint;
    result{i,j,5}=stats;
    m=m+1;
    b1(m)=b(1);
  end
end

bset=result(:,:,1);
a = zeros(n-t,n-t);
b = zeros(n-t,n-t);
flag = n-t;
for j = 1:n-t
    for i = 1:n-t
        if  flag - i  >= 0                                         % only read the upp triangular table
        a(i,j) = bset{i,j}(1);                                     % the slope
        b(i,j) = bset{i,j}(2);                                     % the intercept
        else
            continue;
        end
    end
    flag = flag -1;                                                 % number of rows to  read
end

% slope and intercept
sl = a(a(:) ~=0);
interp = b(b(:) ~=0);

figure
% histogram to a
subplot(1,2,1)
histfit(sl, 20, 'normal')
% mu and sigma
pd_sl = fitdist(sl, 'normal');
title('Regression ratio')

% samples
% R = -betarnd(pd_sl.a, pd_sl.b, [5000,1]);
 R = pd_sl.sigma*randn(5000,1)+pd_sl.mu;

%%
X01=0:(2030-endyear);

x=[year2(:,1)
    X01'];
y=[logprice];

% prediction to 2030
yCalc1 = year2*bbbb;
% yf_time= R*X01+ yCalc1(end);
yf_time= R*X01+ y(end);
%%%%%%
yfmean_time=mean(yf_time);
yfqant_time = quantile(yf_time, [0.05,0.1, 0.25, 0.5,0.75 0.9 0.95]);
yy_time=[yCalc1(1:endyear-startyear)                                               % delete the last year
    yfmean_time'];

figure 
semilogy(startyear:endyear,exp(y),'o')
hold on
semilogy(2030,[1.26 1.39 1.62 2.00 2.4],'o')
hold on
semilogy(startyear:2030,exp(yy_time))
hold on 
semilogy(endyear:2030, exp(yfqant_time), 'b-')
ylim([0.1, 10])
xlim([1970, 2030])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

capacity=wind(:,2);
logcapacity=log(capacity);


year3=[year,ones(n,1)];
[bbbb,bint2,r2,rint2,regress2]=regress(logcapacity,year2);


%%

X01=0:(2030-endyear);
X02=year3(:,1);

x=[X02,ones(n,1)];
ycapacity=[logcapacity];

% prediction to 2050
yCalc1_his = x*bbbb;
yCalc1= bbbb(1,:)*X01+ ycapacity(end);


predictor=cell(n-t,n-t);
predict=cell(n-t,n-t);
for j=1:n-t
for i=1:n-(t-1)-j
  predict{j,i}=logprice(i:n+1-j,:);                         %e
  predictor{j,i}=[logcapacity(i:n+1-j),ones((n+2)-j-i,1)];  
end
end

result={};
b1=zeros((n-(t-1))*(n-t)/2,1);
m=0;
for i=1:n-t
  for j=1:n-(t-1)-i
    [b,bint,r,rint,stats]=regress(predict{i,j},predictor{i,j});
    result{i,j,1}=b;
    result{i,j,2}=bint;
    result{i,j,3}=r;
    result{i,j,4}=rint;
    result{i,j,5}=stats;
    m=m+1;
    b1(m)=b(1);
  end
end

LR=zeros(n-t,n-t);
bset=result(:,:,1);
for i=1:n-t
    for j=1:(n-(t-1))-i
        LR(i,j)=1-2^(bset{i,j}(1));
    end
end
a=LR;
% b=LR(1,:); %  ending year 2006
a(a==0) = [];
Y = prctile(a(:),[5 50 95],1);


%%%%%histogram
figure

histfit(a, 20, 'normal')
title('Learning rate')

% mu and sigma
% pd_sl = fitdist(a', 'normal');

figure

histfit(b1, 20, 'normal')
title('Regression ratio')
pd_sl = fitdist(b1, 'normal');

% samples
R = pd_sl.sigma*randn(5000,1)+pd_sl.mu;


X03=yCalc1(1:2030-endyear+1)-ycapacity(end);

x=yCalc1;  %%log capacity
y=[logprice];

% prediction to 2050

yf= R*X03+y(end);
yfmean=mean(yf);
yfqant = quantile(yf, [0.05,0.1, 0.25, 0.5,0.75 0.9 0.95]);

figure
plot(startyear:endyear,exp(y),'o')
hold on
plot(endyear:2030,exp(yf), 'Color',[0.7,0.7,0.7] )
hold on 
plot(endyear:2030, exp(yfqant), 'r-')
hold on 


xlabel('Capacity')
ylabel('LCOE [2016$/w]')
title('Forecast for solarPV')
legend('historical data','exponential trend projection','Location','southwest')


figure
semilogy(startyear:endyear,exp(y),'o')
hold on
semilogy(endyear:2030,exp(yf), 'Color',[0.7,0.7,0.7] )
hold on 
semilogy(endyear:2030, exp(yfqant), 'r-')
hold on 
y2030 = [min(yf(:, 2030-endyear+1)), max(yf(:, 2030-endyear+1))];
semilogy([2030 2030], exp(y2030), 'k-')
xlabel('Capacity')
ylabel('[2016$/w]')
title('Forecast for solarPV (function of delopyment)')
legend('historical data','exponential trend projection','Location','southwest')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  comparison  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
semilogy(startyear:endyear,exp(y),'o')
hold on
semilogy(startyear:2030,exp(yy_time))
hold on 
semilogy(endyear:2030, exp(yfqant_time), 'b-')
hold on
semilogy(endyear:2030, exp(yfqant), 'g-')

ylim([0.1, 10])
xlim([1970, 2030])

