function wrapper_srtsm(enddate)

load maturities.mat

load data_forwardrates.mat

end_id = find(time ==(200000+round(str2double(enddate))));

rlb = 0.25;

 

%% Lines 35-53

options1  =  optimset('fminunc');

options1  =  optimset(options1 , 'TolFun', 1e-16);

options1  =  optimset(options1 , 'TolX', 1e-16);

options1  =  optimset(options1 , 'Display', 'iter');

options1  =  optimset(options1 , 'Diagnostics', 'off');

options1  =  optimset(options1 , 'LargeScale', 'off');

options1  =  optimset(options1 , 'MaxFunEvals', 10^6) ;

options1  =  optimset(options1 , 'FinDiffType', 'central') ;

options1  =  optimset(options1 , 'MaxIter', 10^6) ;

load parameters.mat

startv = parameters;

epsilon = 1e3;

 

while epsilon > 1e-6

    [parameters,fval] = fminunc(@(parameters)EKF_SRTSM([parameters;0.25],maturities,forwardrates(1:end_id,:),0),startv,options1);

    epsilon = abs(parameters-startv);

    epsilon = max(epsilon);

    startv = parameters;

end

%% Our own code to define maturities

J = (1:120)';

JJ = 1:1:max(J);

JJ = JJ';

%% Lines 10-16 from EKF_SRTSM.m

rhoP = parameters(1:9); rhoP = reshape(rhoP,3,3);

muP = parameters(10:12);

rhoQ1 = parameters(13); rhoQ2 = parameters(14);

sigma = [abs(parameters(15)) 0 0 ; parameters(16) abs(parameters(18)) 0; parameters(17) parameters(19) abs(parameters(20))];

omega = sigma*sigma'; sigma11 = omega(1,1); sigma22 = omega(2,2); sigma33 = omega(3,3); sigma12 = omega(1,2); sigma13 = omega(1,3); sigma23 = omega(2,3);

delta0 = parameters(21);

omegaM = diag(parameters(22)^2*ones(length(J),1));

%% Redundant

rlb = 0.25;

 

%% Lines 26-42 from EKR_SRTSM.m

zt = NaN(length(J),6);

xt = rhoQ1^2; if abs(xt-1) > 1e-5 zt(:,1) = (xt.^J-1)/(xt-1); else zt(:,1) = J; end

xt = rhoQ2^2; if abs(xt-1) > 1e-5 zt(:,2) = (xt.^J-1)/(xt-1); else zt(:,2) = J; end

xt = rhoQ2^2; if abs(xt-1) > 1e-5 zt(:,3) = xt.^J.*(J.^2*xt^2-2*J.^2*xt+J.^2-2*J*xt^2+2*J*xt+xt^2+xt)/(xt-1)^3-xt*(xt+1)/(xt-1)^3; else zt(:,3) = (J-1).*J.*(2*J-1)/6; end

xt = rhoQ1*rhoQ2; if abs(xt-1) > 1e-5 zt(:,4) = (xt.^J-1)/(xt-1); else zt(:,4) = J; end

xt = rhoQ1*rhoQ2; if abs(xt-1) > 1e-5 zt(:,5) = xt/(xt-1)^2 - xt.^J.*(J+xt-J*xt)/(xt-1)^2; else zt(:,5) = J.*(J-1)/2; end

xt = rhoQ2^2; if abs(xt-1) > 1e-5 zt(:,6) = xt/(xt-1)^2 - xt.^J.*(J+xt-J*xt)/(xt-1)^2; else zt(:,6) = J.*(J-1)/2; end

sigmasJ2 = zt*[sigma11; sigma22; (1/rhoQ2^2)*sigma33; 2*sigma12; 2/rhoQ2*sigma13; 2/rhoQ2*sigma23];

sigmasJ = sqrt(sigmasJ2);

bn = [rhoQ1.^JJ rhoQ2.^JJ JJ.*rhoQ2.^(JJ-1)];

Bn = [1 1 0; bn]; Bn = cumsum(Bn,1); Bn = Bn(1:end-1,:);

aJ = NaN(length(J),1);

for j = 1:1:length(J)

    aJ(j) = delta0 - 0.5*Bn(J(j),:)*sigma*sigma'*Bn(J(j),:)'/1200;

end

bJ = [rhoQ1.^J rhoQ2.^J J.*rhoQ2.^(J-1)];

 

% Call EKR_SRTSM.m with index=2 (final arguments) to get Xf

Xf = EKF_SRTSM([parameters;0.25],maturities,forwardrates(1:end_id,:),2);

% Pick the final period's factors

XT = Xf(:,end);

% Create a matrix to hold mean/median forecasts of the factors

XTplusj = zeros(size(Xf,1),120);

% Create a matrix to hold the predicted foreward rates

PredictedForwardRates = zeros(120,120);

% Create a 1 period forecast of forward rates as in lines 45-47 EKR_SRTSM.m

XTplusj(:,1) = muP + rhoP * XT;

musJ = aJ + bJ*XTplusj(:,1);

z1_temp = (musJ-rlb)./sigmasJ;

PredictedForwardRates(:,1) = rlb + (musJ-rlb).*normcdf(z1_temp) + sigmasJ.*normpdf(z1_temp);

% Repeat for 2 to 120 months ahead

for j = 2:120

    % Forecast the factors under the p-measure

    XTplusj(:,j) = muP + rhoP * XTplusj(:,j-1);

    % Price under the Q meaure

    musJ = aJ + bJ*XTplusj(:,j);

    z1_temp = (musJ-rlb)./sigmasJ;

    PredictedForwardRates(:,j) = rlb + (musJ-rlb).*normcdf(z1_temp) + sigmasJ.*normpdf(z1_temp);

end

 

% Converge forward rates to interest rates.

PredictedInterestRates = 0 * PredictedForwardRates;

PredictedInterestRates(1,:) = PredictedForwardRates(1,:);

for j = 2:120

    PredictedInterestRates(j,:) = (PredictedInterestRates(j-1,:)*(j-1) + PredictedForwardRates(j,:)) / j;

end

 

% Save everything

save(['RealTimeEstimates' enddate '.mat'])

end
