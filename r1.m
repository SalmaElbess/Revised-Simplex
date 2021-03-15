% implementation of phase 1 revised simplex based on intoroduction to
%operations research by hiller
clc;clear;
%% problem initialization
%Wyndor Glass Co. problem from intoroduction to operations research by
%hiller
c_initial = [3,5,0,0,0];
b = [4;12;18]; %b
A = [1,0;0,2;3,2]; % A matrix
x = [1;2;3;4;5]; %variables
m = length(b); %number of constrains
A_I = [A eye(m)]; % A matrix for all variables including the slacke ones
B_inv_initial = eye(m);
xs = [3,4,5]; %slack_variables
xb = xs; %inital basic variables

xb_values = b;
B_inv = B_inv_initial;
c = - c_initial;

n_iterations = 0;
while(true)
    n_iterations = n_iterations +1;
    %% choose entring variable
    [value,xk] = min(c);%entering basic variable
    aik = A_I(:,xk); %coefficients of xk
    %% choose leaving variable
    [val,r] = get_min_postivive(xb_values./aik); % r = number of equation containing the leaving basic variable
    xb(r)= xk; %switch between entering and leaving variables in basic variables list
    %% update B_inv
    ark = aik(r);
    E = eye(m);
    eta = -1*aik./ark;
    eta(r) = 1/ark;
    E(:,r) = eta;
    B_inv = E*B_inv;
    %% update basic variables values
    xb_values = B_inv*b;
    %% optimality test
    cb = c_initial(xb); %weights of basic variables in f equation
    c =  cb*B_inv*A_I - c_initial;
    if check_optimality(c)
        break;
    end
end
objective_f_value = cb*B_inv*b;
%% matlab built in function 
x = linprog(-[3,5],A,b);
objective_built_in = [3,5]*x;
function optim = check_optimality(c)
%check that all the variables are non negative
if any(c <0)
    optim = false;
else
    optim = true;
end
end
function [val,index] = get_min_postivive(c)
% get the min positive value and index in an array
val = inf;
index = 0;
for i = 1:length(c)
    if c(i) < val && c(i) > 0
        val = c(i);
        index = i;
    end
end
end