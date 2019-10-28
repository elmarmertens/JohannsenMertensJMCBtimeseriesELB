[~, ~, datacell] = xlsread('feds200628.xls');
datacell = datacell(11:end,:);
try
    date = datenum(datacell(:,1),'mm/dd/yyyy');
catch
    date = datenum(datacell(:,1),'yyyy-mm-dd');
end
Y = year(date); M = month(date); D = day(date);
% Pick the data from 1990/01, note that they are in reverse
% chronological order.
id_end = find(100*Y+M == 199001,1,'last');
Y = Y(1:id_end); M = M(1:id_end); D = D(1:id_end); data = cell2mat(datacell(1:id_end,end-5:end));
GSWdata = [Y M D data];GSWdata = flipud(GSWdata);
% Find the end of month observation
[b,m,n] = unique(100*GSWdata(:,1)+GSWdata(:,2),'last');
GSWdata = GSWdata(m,:); 
% Get rid of the dates without any parameter values
flag = isnan(GSWdata(:,4:end)); flag = sum(flag,2); id = find(flag == 6);
GSWdata(id,:) = [];
% Determine whether the last observation is the end of month 
if GSWdata(end,3) < day(lbusdate(GSWdata(end,1),GSWdata(end,2)))
    GSWdata(end,:) = [];
end
time = GSWdata(:,1)*100 + GSWdata(:,2);

%Calculate yields from Svensson formula
beta0 = GSWdata(:,4)'; beta1 = GSWdata(:,5)'; beta2 = GSWdata(:,6)'; beta3 = GSWdata(:,7)'; tau1 = GSWdata(:,8)'; tau2 = GSWdata(:,9)';
n = (1:120)'/12;
beta0 = repmat(beta0,length(n),1); beta1 = repmat(beta1,length(n),1); beta2 = repmat(beta2,length(n),1); beta3 = repmat(beta3,length(n),1); tau1 = repmat(tau1,length(n),1); tau2 = repmat(tau2,length(n),1);
n = repmat(n,1,size(GSWdata,1));
yields = beta0 + beta1.*(1-exp(-n./tau1))./n.*tau1 + beta2.*((1-exp(-n./tau1))./n.*tau1-exp(-n./tau1)) + beta3.*((1-exp(-n./tau2))./n.*tau2-exp(-n./tau2));
yields = [yields NaN(size(yields,1),24)];

save('xlsdata.mat','time','yields');
