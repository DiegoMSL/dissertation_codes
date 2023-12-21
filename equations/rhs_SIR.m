% -----------------------------------------------------------------
% EPIDEMIC - Epidemiology Educational Code
% www.EpidemicCode.org
% -----------------------------------------------------------------
% This function defines the system of ODEs for the
% SIR epidemic model.
%
% The dynamic state coordinates are:
%
%  S = susceptibles        (number of individuals)
%  I = infected            (number of individuals)
%  R = recovered           (number of individuals)
%
% The epidemic model parameters are:
%
%   N     = population size   (number of individuals)
%   beta  = transmission rate (days^-1)
%   gamma = recovery rate     (days^-1)
% -----------------------------------------------------------------
% programmers: Eber Dantas
%              Americo Cunha
%
% last update: May 19, 2020
% -----------------------------------------------------------------

% -----------------------------------------------------------------
function dydt = rhs_SIR(t,y,param)

% model parameters: param = [N beta gamma]
N     = param(1);  % population size   (number of individuals)
beta  = param(2);  % transmission rate (days^-1)
gamma = param(3);  % recovery rate     (days^-1)

% SIR dynamic model:
% 
%      y = [S I R]             is the state vector
%   dydt = [dSdt dIdt dRdt] is the evolution law
% 
% dSdt - rate of susceptible         (number of individuals/days)
% dIdt - rate of infected            (number of individuals/days)
% dRdt - rate of recovered           (number of individuals/days)

[S I R] = deal(y(1),y(2),y(3));

dSdt = - beta*S.*(I/N);
dIdt = beta*S.*(I/N) - gamma*I;
dRdt = gamma*I;
dydt = [dSdt; dIdt; dRdt];

end
% -----------------------------------------------------------------