% BMED 3310 Computational Project
% Part 1 Eva Beresin Spring 2026
clear; close all; clc;

% Load data under my name key (175210)
projectdata = ProjectDataGenerator(175210);

% Create a dimensional matrix
% Rows: mass (M), length (L), time (T) Columns 1-3 (basis): length 1, mass diffusivity 1, pressure gradient 1 Columns 4-9 (non-basis): length 2, mass resistance 1, mass 1, velocity 1, power density 1, mass resistance 2 Thermal conductivity 1 is removed because the basis parameters don't cover temperature (Kelvin)

er_matrix = [
    0, 0, 1, 0, 0, 1, 0, 1, 0; % mass exponents
    1, 2, -2, 1, -3, 0, 1, -1, -3; % length exponents
    0, -1, -2, 0, 1, 0, -1, -3, 1  % time exponents
];

% Rearrange dimensional matrix and use rref to obtain the pi exponent table
eva_data = rref(er_matrix);
disp('Pi Exponent Table:');
disp(eva_data);

% Calculate the pi groups from observation data
% Extract the 200x10 measurements table
measurements = projectdata{2,2};

% Extract the basis parameters by column index
b1 = measurements{:, 1}; % length 1
b2 = measurements{:, 9}; % mass diffusivity 1
b3 = measurements{:, 7}; % pressure gradient 1

% Extract the non-basis parameters by column index
p1 = measurements{:, 2}; % length 2 (dependent)
p2 = measurements{:, 3}; % mass resistance 1
p3 = measurements{:, 4}; % mass 1
p4 = measurements{:, 5}; % velocity 1
p5 = measurements{:, 6}; % power density 1
p6 = measurements{:, 8}; % mass resistance 2

% Form the pi groups using the pi exponent table (rr_matrix)
% Formula to use: pi = (b1^a .* b2^b .* b3^c) ./ p
% The exponents a, b, and c come from columns 4 through 9 of the rr matrix

pi_1 = (b1.^1 .* b2.^0 .* b3.^0) ./ p1;
pi_2 = (b1.^-1 .* b2.^-1 .* b3.^0) ./ p2;
pi_3 = (b1.^6 .* b2.^-2 .* b3.^1) ./ p3;
pi_4 = (b1.^-1 .* b2.^1 .* b3.^0) ./ p4;
pi_5 = (b1.^-1 .* b2.^1 .* b3.^1) ./ p5;
pi_6 = (b1.^-1 .* b2.^-1 .* b3.^0) ./ p6;

% Combine into one table of pi group values (200 rows x 6 columns)
pi_data = [pi_1, pi_2, pi_3, pi_4, pi_5, pi_6];

% Output to the command window in scientific notation
format short e;
disp('Table of Pi Group Values:');
disp(pi_data);


% BMED 3310 Computational Project
% Part 2 Eva Beresin Spring 2026

% Scale the pi group data
% normalize to max (element-wise division so each column is scaled by its own max)
pi_scaled = pi_data ./ max(pi_data);

% Isolate the dependent and independent variables
% dependent parameter is Length 2 --> pi_1 (column 1)
dep_y = pi_scaled(:, 1);

% The remaining columns (2-6) are independent predictor variables
indep_x = pi_scaled(:, 2:end);

% Multi-linear fit using fitlm
disp('--- Multi-linear Model (fitlm) ---');
multimodel = fitlm(indep_x, dep_y, 'linear');
disp(multimodel);

% Multi-linear fit using fitnlm
% model: y = b(1) + b(2)*x1 + b(3)*x2 + b(4)*x3 + b(5)*x4 + b(6)*x5
modelfunc = @(b,x) b(1) + x*b(2:end);

% Initial guess --> coefficients
initialguess = multimodel.Coefficients.Estimate;

disp('--- Multi-linear Model (fitnlm) ---');
% do fit
multimodel2 = fitnlm(indep_x, dep_y, modelfunc, initialguess);
disp(multimodel2);

% Exponential fit using fitnlm
% model: y = b(1)*exp(b(2)*x1 + b(3)*x2 + ...)
modelexp = @(b,x) b(1) * exp(x*b(2:end));

% initial guess for exponential
% linear fit of independent data vs. log(dependent data)
linexp = fitlm(indep_x, log(dep_y), 'linear');
nextguess = linexp.Coefficients.Estimate;

% update scaling coefficient to be the exponential of the estimate
nextguess(1) = exp(nextguess(1));

disp('--- Exponential Model (fitnlm) ---');
% do fit
expmodel = fitnlm(indep_x, dep_y, modelexp, nextguess);
disp(expmodel);

% Power fit using fitnlm
% model: y = b(1)*x1^b(2)*x2^b(3)*x3^b(4)*x4^b(5)*x5^b(6)
powermodel = @(b,x) b(1) * (x(:,1).^b(2)) .* (x(:,2).^b(3)) .* (x(:,3).^b(4)) .* (x(:,4).^b(5)) .* (x(:,5).^b(6));

% Initial guess for power: linear fit of log(independent data) vs. log(dependent data)
linfitpower = fitlm(log(indep_x), log(dep_y), 'linear');
powerguess = linfitpower.Coefficients.Estimate;

% Update the scaling coefficient to be the exponential of the estimate
powerguess(1) = exp(powerguess(1));

disp('--- Power Model (fitnlm) ---');
% Do fit
finalmodel = fitnlm(indep_x, dep_y, powermodel, powerguess);
disp(finalmodel);


% BMED 3310 Computational Project
% Part 3 Eva Beresin Spring 2026
% get the max values of the original unscaled pi groups from part 1
% pi_data has 6 columns. column 1 is dependent, columns 2-6 are independent

ymax = max(pi_data(:, 1));
xmax = max(pi_data(:, 2:end)); % 1x5 array of max independent values

% unscaling coefficients
% unscaling multi-linear coefficients
unscaledml = zeros(6, 1);
unscaledml(1) = ymax * multimodel2.Coefficients.Estimate(1);
unscaledml(2:end) = (ymax ./ xmax') .* multimodel2.Coefficients.Estimate(2:end);

disp('Unscaled Multi-Linear Coefficients');
disp(unscaledml);

% unscaling exponential coefficients
unscaledexp = zeros(6, 1);
unscaledexp(1) = ymax * expmodel.Coefficients.Estimate(1);
unscaledexp(2:end) = expmodel.Coefficients.Estimate(2:end) ./ xmax';

disp('Unscaled Exponential Coefficients');
disp(unscaledexp);

% unscaling power coefficients
powerunscaled = zeros(6, 1);
bpower = finalmodel.Coefficients.Estimate;

% new intercept --> ymax * b0 / (xmax1^b1 * xmax2^b2 * ...)
powerunscaled(1) = ymax * bpower(1) / prod(xmax'.^bpower(2:end));
powerunscaled(2:end) = bpower(2:end); % exponents don't change

disp('Unscaled Power Coefficients');
disp(powerunscaled);

% choosing the best fit
% get AIC values
multi_aic = multimodel2.ModelCriterion.AIC;
exp_aic = expmodel.ModelCriterion.AIC;
pow_aic = finalmodel.ModelCriterion.AIC;

fprintf('\n Choosing Best Fit \n');
fprintf('The AIC value of the multilinear model is %f\n', multi_aic);
fprintf('The AIC value of the exponential model is %f\n', exp_aic);
fprintf('The AIC value of the power model is %f\n', pow_aic);

% the best model has the lowest/most negative AIC
[~, bestmodel] = min([multi_aic, exp_aic, pow_aic]);
models = {'Multi-Linear', 'Exponential', 'Power'};
fprintf('The best model is %s because it has the lowest AIC value.\n', models{bestmodel});

% get critical pi groups
% initial significance level alpha = 0.01
% bonferroni correction alpha_adjusted = alpha / n (where n=5 independent variables)
adjustedalpha = 0.01 / 5;

fprintf('\n Identify Critical Pi Groups \n');
fprintf('Critical pi groups were identified based on p-value < %f\n', adjustedalpha);

% get out p-values from the best exponential model
chosenpvalues = expmodel.Coefficients.pValue(2:end);

% get a list of all the non-basis parameters for cols 2-6
parameters = {'Mass Resistance 1', 'Mass 1', 'Velocity 1', 'Power Density 1', 'Mass Resistance 2'};

fprintf('The non-basis parameters for the critical pi groups are: \n');
for i = 1:5
    if chosenpvalues(i) < adjustedalpha
        fprintf('- %s\n', parameters{i});
    end
end
