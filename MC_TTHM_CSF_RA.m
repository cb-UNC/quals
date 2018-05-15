% This code is written by me, Christopher Bowers, for my qualifying exam on
% 05/14/18 . This is to answer question 5 of the exam.
%
% The main function below carries out a Monte Carlo simulation which takes
% inputs from a multivariate distribution in order to calculate the
% total trihalomethane (TTHM) present in a produced water with a source water
% bromide concentration between 20 and 40 ug/L. It also calculates the
% lifetime probability of cancer for consumers of this water due to TTHM
% exposure. Here, TTHM refers specifically to TTHM4.
%
% The indexing scheme for this system refers to different THMs by:
%           1   -   Chloroform
%           2   -   BDCM
%           3   -   DBCM
%           4   -   Bromoform
%
% Function is Monte-Carlo Total TriHaloMethane Cancer Slope Factor Risk 
% Assessment (MC_TTHM_CSF_RA) 

function [mean_CR,var_CR,mu_CR,sigma_CR,R,CR_PLOT,P_CR_PLOT,CR_TARGET_I,P_TARGET_CR,TTHM4_i_80,P_80,P_TTHM4_PLOT,TTHM4_PLOT,MU_MC,SIGMA_MC,MEAN_MC,p_MC] = MC_TTHM_CSF_RA
% Model inputs - lognormal Means
a1 = 9.77; 
a2 = 8.36; 
a3 = 5.79; 
a4 = 4.98;
MU = [a1 a2 a3 a4];

% Model inputs - lognormal standard deviations
b = zeros(1,4);
b(1) = 1.01;
b(2) = 2.35;
b(3) = 2.98;
b(4) = 2.03;

% Model inputs - correlation factors
p = zeros(4,4);
p(1,2) = 0.613;     p(1,1) = 1;     p(2,1) = p(1,2);
p(1,3) = 0.475;     p(2,2) = 1;     p(3,1) = p(1,3);
p(1,4) = 0.443;     p(3,3) = 1;     p(4,1) = p(1,4);
p(2,3) = 0.501;     p(4,4) = 1;     p(3,2) = p(2,3);
p(2,4) = 0.335;                     p(4,2) = p(2,4);
p(3,4) = 0.231;                     p(4,3) = p(3,4);

% Generate covariance matrix from inputs
COV = zeros(4,4);
for i=1:4
    for j=1:4
        COV(i,j) = p(i,j)*b(i)*b(j);
    end
end

% Randomly sample from the multivariate normal distribution N times
N = 100000;
R = mvnrnd(MU,COV,N);

% Check the normal (i.e. lognormal) mean and standard deviation
MU_MC = mean(R);
SIGMA_MC = var(R);
for i=1:4, SIGMA_MC(i) = SIGMA_MC(i)^(1/2); end

% Caclulate the sample mean for each THM
MEAN_MC = zeros(1,4);
for i=1:4, MEAN_MC(i) = exp(MU_MC(i) + (SIGMA_MC(i))^2/2)*(1e-3); end

% Calculate the correlation coefficients for comparison
p_MC = zeros(4,4);
for i=1:4
    for j=1:4
        cov_ij = 0;
        for k=1:N
            cov_ij = cov_ij + (R(k,i)-MU_MC(i))*(R(k,j)-MU_MC(j))/N;
        end
        p_MC(i,j) = cov_ij/(SIGMA_MC(i)*SIGMA_MC(j));
    end
end

% Generate the empirical CDF for TTHM4
TTHM4 = arrayfun(@(x) exp(x), R);       % First, convert from ln(Y) to X
TTHM4 = TTHM4 .* 10^(-3);               % Convert to ug/L
TTHM4 = sum(TTHM4,2);                   % Sum the THM concentrations to get TTHM4
TTHM4 = sort(TTHM4);                    % Sort to get the empirical CDF

% Generate probability space
P = linspace(1,N,N);
P = P./N;

% Generate plottable arrays which are of a manageable size for export, and 
% calculate probability that the samples are below 80 ug/L while you're at
% it.
TTHM4_MAX = 300;
TTHM4_80 = 80;
i_80 = 0;
for i=1:N
    if (TTHM4(i) > TTHM4_MAX)
        i_END = i;
        break
    end
    if (TTHM4(i) > TTHM4_80) && (i_80 == 0)
        i_80 = i;
    end
end

% Calculate the probability that the concentration is less than 80 ug/L
P_80 = P(i_80);
TTHM4_i_80 = TTHM4(i_80);

% Cut the plotted points down by a factor of 100 for simplicity
i_PLOT = 1:100:i_END;
P_TTHM4_PLOT = P(i_PLOT);
TTHM4_PLOT = TTHM4(i_PLOT);

% Plot the ECDF of TTHM4 leaving the treatment plants
figure
plot(TTHM4_PLOT,P_TTHM4_PLOT)



% Calculate the TTHM4 cancer risk
CR = arrayfun(@(x) exp(x), R);       % First, convert from ln(Y) to X
CR = CR .* 10^(-6);                  % Convert to mg/L
h = 2/70;   CR = CR .* h;            % Calculate LADDs

% Cancer risk from Chloroform
h = 6.1e-3; CR(:,1) = CR(:,1) .* h;

% Cancer risk from Bromodichloromethane
h = 6.2e-2; CR(:,2) = CR(:,2) .* h;

% Cancer risk from Dibromochloromethane
h = 8.4e-2; CR(:,3) = CR(:,3) .* h;

% Cancer risk from Bromoform
h = 7.9e-3; CR(:,4) = CR(:,4) .* h;

% Total accumulated cancer risk
CR = sum(CR,2);

% Generate ECDF
CR = sort(CR);
CR_MAX = 2.5e-4;
CR_TARGET = 10^(-5);
i_TARGET = 0;
for i=1:N
    if (CR(i) > CR_MAX)
        i_END = i;
        break
    end
    if (CR(i) > CR_TARGET) && (i_TARGET == 0)
        i_TARGET = i;
    end
end

P_TARGET_CR = P(i_TARGET);
CR_TARGET_I = CR(i_TARGET);

% Cut the plotted points down by a factor of 100 for simplicity
P = P.';
i_PLOT = 1:100:i_END;
P_CR_PLOT = P(i_PLOT);
CR_PLOT = CR(i_PLOT);

% Plot the ECDF of the cancer risk associated with TTHM4
figure
plot(CR_PLOT,P_CR_PLOT)

%

% Now we will calculate the mean and standard deviation of the lifetime
% cancer risk associated with DBP exposure. We do this by fitting a
% log-normal CDF to our empirical CDF
fun = @(x,xdata) 0.5*(1 + erf((log(xdata)-x(1))/(x(2)*(2)^(1/2)) ) );
x0 = [9,0.25]; x0(1) = log(mean(CR))-x0(2)/2;
x = lsqcurvefit(fun,x0,CR,P);
mu_CR = x(1); sigma_CR = x(2);

% Plot the fitted values with the ECDF to compare
figure
plot(CR_PLOT,P_CR_PLOT,'ko',CR_PLOT,fun(x,CR_PLOT),'b-');

mean_CR = mean(CR); var_CR = var(CR);

end