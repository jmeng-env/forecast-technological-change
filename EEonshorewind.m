%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% onshore wind (investment cost)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c= [  0 0 1;
      0 1 0;
      1 0 0;
      1 1 0;
      1 0 1;
      0 1 1;
      0.4 0.4 0.4;
      0.7 0.1 0.2;
      0.9 0.5 0.2;
      0.1 0.9 0.9];
  
num2014 = xlsread('Wind.xlsx','EE','A6:E6');
num_new2014 = num2014*1.023967134;  %2014 price to 2016
num_new2014(:,1) = num2014(:,4);
num_new2014(:,2) = num2014(:,2);
num_new2014(:,3) = num2014(:,1);
num_new2014(:,4) = num2014(:,3);
num_new2014(:,5) = num2014(:,5);

num2030 = xlsread('Wind.xlsx','EE','A10:E10');
num_new2030 = num2030*1.023967134;  %2010 price to 2016
num_new2030(:,1) = num2030(:,4);
num_new2030(:,2) = num2030(:,2);
num_new2030(:,3) = num2030(:,1);
num_new2030(:,4) = num2030(:,3);
num_new2030(:,5) = num2030(:,5);

 t=[1 3];


for i=1:5
    a(i)=log(num_new2030(i)/num_new2014(i))/(16-1);
    b(i)=log((num_new2014(i)^16)/num_new2030(i))/(16-1);
end

for i=1:2
    for j=1:5
         EE(i,j)=exp(a(j)*t(i)+b(j));
    end
end
num= xlsread('D:\OneDrive - University Of Cambridge\Cambridge\INNOPATH\compare\regressiondata.xlsx','2015','U6:Y6');   % LR
for i=1
  num_new(i,2) = num(i,1);
  num_new(i,3) = num(i,3);
  num_new(i,4) = num(i,5);
  num_new(i,1) = num(i,1)*num(i,1)/num(i,3);
  num_new(i,5) = num(i,5)*num(i,5)/num(i,3);
end

% check the skewness
sk = (EE(:,1) - EE(:,2))./(EE(:,3) - EE(:,2));
EE(10,:)=num_new;

% using the smoothingsplines spline
yData = [0; 0.1; 0.5; 0.9; 1];

figure
xData = EE';

N = 50000;

fx = zeros(N-1, size(EE,1));
fitresult = zeros(N, size(EE,1));


for i =[2]
%     fit the data betwn the min and max
    x_range(i,:) = linspace(min(min(EE(i,:))), max(max(EE(i,:))), 50000);
    idx1 = find(x_range(i,:) >= min(xData(:,i)) & x_range(i,:)  <= max(xData(:,i)));
    fitresult(idx1,i) =  interp1(xData(:,i), yData, x_range(i,idx1), 'pchip' );
    
%     pad the data on the right (Note: Not nececssary to pad the left side since it is already 0)
    idx2= find(x_range(i,:) > max(xData(:,i)) );
    fitresult(idx2,i) = 1;
    
   
%     Plot fit with data.
    hold on
    h = plot(x_range(i,:), fitresult(:,i), 'Color', c(i,:));
    h.LineWidth = 1;
    h2 = plot(xData(:,i), yData, 'rx');
    h2.LineWidth = 2;
    
    h = (x_range(i,2)-x_range(i,1));
    fx(:,i) = diff(fitresult(:,i))./h; % approximation
%     pause
end

xlim([0, 1])    

xlabel('LCOE for wind (2016US$w)')
ylabel('Cumulative Probability')
grid on
hold off
legend('2017','Location','northeast')



xlabel(' investment cost for US utility solar PV (2016US$w)')
ylabel('Probability Density')
grid on
axis tight
hold off
legend('2015','Location','northeast')

% refit a distribution to the pdf

for i = 2
      h5 = plot(x_range(i,2:end), fx(:,i), '.-');
%     fit a gamma distribution to the empircal pdf
    modelFun3 =  @(p,x) (p(2).^p(1)) .*  x.^(p(1)-1) .* exp(-p(2).*x) ./gamma(p(1));  % gamma

    startingVals = [7 1];
    coefEsts3 = nlinfit(x_range(i,2:end), fx(:,i)', modelFun3, startingVals);
    xgrid = linspace(0, max(x_range(i,:)),5000);
    line(xgrid, modelFun3(coefEsts3, xgrid), 'Color', c(1,:));
    hold on
end
q3 = gamcdf(1.38,coefEsts3(1),1/coefEsts3(2));


xlabel('investment cost for onshore wind (2016US$W)')
ylabel('Probability Density')
grid on
axis tight
legend('2015','Location','northeast')
y=0:3;
x1=1.38*ones(1,4);  %global 2016
plot(x1,y, 'Color',c(10,:));
