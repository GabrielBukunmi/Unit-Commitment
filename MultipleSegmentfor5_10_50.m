clear; clc; close all
syms p1 p2 p3

h1=510+7.2*p1+0.00142*p1.^2;
h2=310+7.85*p2+0.00194*p2.^2;
h3=78+7.97*p3+0.00482*p3.^2;

%Fuel cost for each generator
fuel_cost1=1.1;
fuel_cost2=1.0;
fuel_cost3=1.2;

%Get the cost function
 F1=fuel_cost1*h1;
 F2=fuel_cost2*h2;
 F3=fuel_cost3*h3;
% generator 1
p1_min = 150; 
p1_max = 600;
% generator 2
p2_min = 100;
p2_max = 400;
% generator 3
p3_min = 50;
p3_max = 200;


unit1=0;

segment_division = 50; %Ask user to enter the number of segments
disp(['Using ' num2str(segment_division) ' segments']);
results = cell(segment_division, 1); % Cell array to store results

for seg = 1:segment_division
    step1 = (p1_max - p1_min) / segment_division;
    step2 = (p2_max - p2_min) / segment_division;
    step3 = (p3_max - p3_min) / segment_division;

    % Initialize values
    p1_values = [];
    p2_values = [];
    p3_values = [];

    current_p1 = p1_min;
    current_p2 = p2_min;
    current_p3 = p3_min;

    while current_p1 <= p1_max && current_p2 <= p2_max && current_p3 <= p3_max
        p1_values = [p1_values, current_p1];
        p2_values = [p2_values, current_p2];
        p3_values = [p3_values, current_p3];

        current_p1 = current_p1 + step1;
        current_p2 = current_p2 + step2;
        current_p3 = current_p3 + step3;
    end

    % Store the generated values for this segment
    results{seg}.p1 = p1_values;
    results{seg}.p2 = p2_values;
    results{seg}.p3 = p3_values;
end


% Initialize arrays to store cost function values
F_1 = [];
F_2 = [];
F_3 = [];

for i = 1:length(p1_values)
    % Calculate cost function values for each generator at p1_values(i)
    F_P_min1 = double(subs(F1, p1, p1_values(i)));

    % Store the calculated values
    F_1 = [F_1, F_P_min1];

    % Calculate cost function values for each generator at p1_values(i)
    F_P_min2 = double(subs(F2, p2, p2_values(i)));

    % Store the calculated values
    F_2 = [F_2, F_P_min2];

    % Calculate cost function values for each generator at p1_values(i)
    F_P_min3 = double(subs(F3, p3, p3_values(i)));

    % Store the calculated values
    F_3 = [F_3, F_P_min3];
end

% Initialize arrays to store slopes
slope_1 = [];
slope_2 = [];
slope_3 = [];

for i = 1:(length(p1_values) - 1)
    % Calculate slopes for each generator at p1_values(i)
    slope_1 = [slope_1, (F_1(i+1) - F_1(i)) / (p1_values(i+1) - p1_values(i))];
    slope_2 = [slope_2, (F_2(i+1) - F_2(i)) / (p2_values(i+1) - p2_values(i))];
    slope_3 = [slope_3, (F_3(i+1) - F_3(i)) / (p3_values(i+1) - p3_values(i))];
end


b1 = [];
b2 = [];
b3 = [];

for i = 1:length(p1_values)
    % Calculate b for each generator at the minimum point of each segment
    if i < length(p1_values)
        b1 = [b1, F_1(i) - slope_1(i) * p1_values(i)];
        b2 = [b2, F_2(i) - slope_2(i) * p2_values(i)];
        b3 = [b3, F_3(i) - slope_3(i) * p3_values(i)];
    end
end

f1 = [];
f2 = [];
f3 = [];

for i = 1:length(slope_1) % Assuming F1, F2, F3 are defined arrays
    f1 = [f1, slope_1(i)];
    f2 = [f2, slope_2(i)];
    f3 = [f3, slope_3(i)];
end


p1 = [];
p2 = [];
p3 = [];
for i = 2:length(p1_values)
    p1 = [p1, p1_values(i)];
    p2 = [p2, p2_values(i)];
    p3 = [p3, p3_values(i)];

end


load = 550; %Ask user to enter the value of the load in MW
%% 

Gen_1=0;


Gen_2=0;

Gen_3=0;
disp('0 means off and 1 means on');

Gen_1= input('Enter if generator 1 is On or Off: ');
Gen_2= input('Enter if generator 2 is On or Off: ');
Gen_3= input('Enter if generator 3 is On or Off: ');
%%  Gen 1,2 and 3 OFF

if(Gen_1==0 && Gen_2 ==0 && Gen_3==0)
    disp('infeasible solution')
end
%% Gen 1,2 OFF and 3 ON

if(Gen_1==0 && Gen_2 ==0 && Gen_3==1)
    disp('infeasible solution')
end

%% Gen 1 OFF, 2 ON and 3 OFF

if(Gen_1==0 && Gen_2 ==1 && Gen_3==0)
    disp('infeasible solution')
end
%% Gen 1 OFF, 2 ON and 3 ON

if(Gen_1==0 && Gen_2 ==1 && Gen_3==1)
     p_min_total = p2_min + p3_min;
    f = [slope_2 slope_3];
Aeq = [ones(size(f))];
beq = 550-p_min_total;
ub = [step2*ones(1,segment_division) step3*ones(1,segment_division)];
lb = zeros((size(f)));
[sol,] = linprog(f, [], [], Aeq, beq, lb, ub);
P2=0; % get final value of three generations 
P3=0;
for i=1:segment_division

    P2=P2+sol(i,1);
        P3=P3+sol(i+segment_division,1);
       
end
P1=0;
P2=P2+p2_min;
P3=P3+p3_min;
 %Calculate the minimum cost
 syms p1 p2 p3
    F1 = 0;
    F2 = 310 + 7.85 * p2 + 0.00194 *  p2.^2; %cost function2
     F3 = 93.6 + 9.564* p3 + 0.005784 *  p3^2;  % cost function3

    Minimumcost = double(subs(F2, p2, P2)) + double(subs(F3, p3, P3));
    %display the final results
   MaxGen=p2_max +  p3_max;
    disp(['Max Gen: ' num2str(MaxGen) ' MW']);
      disp(['Min Gen: ' num2str(p_min_total) ' MW']);
      disp(['P1: ' num2str(P1)]);
   disp(['Optimal P2: ' num2str(P2) ' MW']);
    disp(['Optimal P3: ' num2str(P3) ' MW']);
        disp(['F1: ' num2str(F1)]);
    disp(['F2: ' num2str(double(subs(F2, p2, P2))) ' per hour']);
  disp(['F3: ' num2str(double(subs(F3, p3, P3))) ' per hour']);
    disp(['Total Minimum Cost: $' num2str(Minimumcost) ' per hour']);
end
%% Gen 1 ON, 2 OFF and 3 OFF

if(Gen_1==1 && Gen_2 ==0 && Gen_3==0)
    p_min_total = p1_min;
    f = [slope_1];
Aeq = [ones(size(f))];
beq = 550-p_min_total;
ub = [step1*ones(1,segment_division)];
lb = zeros((size(f)));
[sol,] = linprog(f, [], [], Aeq, beq, lb, ub);
P1=0;
P2=0; % get final value of three generations 
P3=0;
for i=1:segment_division

    P1=P1+sol(i,1);
       
end
P1=P1+p1_min;
P2=0;
P3=0;
 %Calculate the minimum cost
 syms p1 p2 p3
    F1 = 561 + 7.92 * p1 + 0.001562 * p1.^2;
    F2 = 0; 
     F3 =0;
    Minimumcost = double(subs(F1, p1, P1));
    %display the final results
   MaxGen=p1_max;
    disp(['Max Gen: ' num2str(MaxGen) ' MW']);
      disp(['Min Gen: ' num2str(p_min_total) ' MW']);
      disp(['Optimal P1: ' num2str(P1) 'MW']);
   disp(['P2: ' num2str(P2)]);
    disp(['P3: ' num2str(P3) ]);
        disp(['F1: ' num2str(double(subs(F1, p1, P1))) ' per hour']);
    disp(['F2: ' num2str(F2) ]);
  disp(['F3: ' num2str(F3) ]);
    disp(['Total Minimum Cost: $' num2str(Minimumcost) ' per hour']);
end
%% Gen 1 ON, 2 OFF and 3 ON

if(Gen_1==1 && Gen_2 ==0 && Gen_3==1)
    p_min_total = p1_min+p3_min;
    f = [slope_1 slope_3];
Aeq = [ones(size(f))];
beq = 550-p_min_total;
ub = [step1*ones(1,segment_division) step3*ones(1,segment_division)];
lb = zeros((size(f)));
[sol,] = linprog(f, [], [], Aeq, beq, lb, ub);
P1=0;
P2=0; % get final value of three generations 
P3=0;
for i=1:segment_division

    P1=P1+sol(i,1);
        P3=P3+sol(i+segment_division,1);
end
P1=P1+p1_min;
P2=0;
P3=P3+p3_min;
 %Calculate the minimum cost
 syms p1 p2 p3
    F1 = 561 + 7.92 * p1 + 0.001562 * p1.^2;
    F2 = 0; 
     F3 =93.6 + 9.564* p3 + 0.005784 *  p3^2;
    Minimumcost = double(subs(F1, p1, P1))+ double(subs(F3, p3, P3));
    %display the final results
   MaxGen=p1_max + p3_max;
    disp(['Max Gen: ' num2str(MaxGen) ' MW']);
      disp(['Min Gen: ' num2str(p_min_total) ' MW']);
      disp(['Optimal P1: ' num2str(P1) 'MW']);
   disp(['P2: ' num2str(P2)]);
    disp(['Optimal P3: ' num2str(P3) ' MW']);
        disp(['F1: ' num2str(double(subs(F1, p1, P1))) ' per hour']);
    disp(['F2: ' num2str(F2) ]);
     disp(['F3: ' num2str(double(subs(F3, p3, P3))) ' per hour']);
    disp(['Total Minimum Cost: $' num2str(Minimumcost) ' per hour']);
end
%% Gen 1 and 2 ON, Gen 3 OFF


if(Gen_1==1 && Gen_2 ==1 && Gen_3==0)
    p_min_total = p1_min + p2_min;
    f = [slope_1 slope_2];
Aeq = [ones(size(f))];
beq = 550-p_min_total;
ub = [step1*ones(1,segment_division) step2*ones(1,segment_division)];
lb = zeros((size(f)));
[sol,] = linprog(f, [], [], Aeq, beq, lb, ub);
P1=0; % get final value of three generations 
P2=0;
for i=1:segment_division
    P1=P1+sol(i,1);
   P2=P2+sol(i+segment_division,1);
end
P1=P1+p1_min;
P2=P2+p2_min;
P3=0;
 %Calculate the minimum cost
 syms p1 p2
    F1 = 561 + 7.92 * p1 + 0.001562 * p1.^2; %cost funciton1
    F2 = 310 + 7.85 * p2 + 0.00194 *  p2.^2; %cost function2
    F3=0;
    Minimumcost = double(subs(F1, p1, P1)) + double(subs(F2, p2, P2));
    %display the final results
   MaxGen=p1_max +  p2_max;
    disp(['Max Gen: ' num2str(MaxGen) ' MW']);
      disp(['Min Gen: ' num2str(p_min_total) ' MW']);
   disp(['Optimal P1: ' num2str(P1) ' MW']);
    disp(['Optimal P2: ' num2str(P2) ' MW']);
    disp(['P3: ' num2str(P3)]);
    disp(['F1: ' num2str(double(subs(F1, p1, P1))) ' per hour']);
    disp(['F2: ' num2str(double(subs(F2, p2, P2))) ' per hour']);
    disp(['F3: ' num2str(F3)]);
    disp(['Total Minimum Cost: $' num2str(Minimumcost) ' per hour']);
end

%% 

% Gen 1,2 and 3 ON
if(Gen_1==1 && Gen_2 ==1 && Gen_3==1)
    p_min_total = p1_min + p2_min + p3_min;
    %formulate the linprog inputs
    f = [slope_1 slope_2 slope_3];
    Aeq = [ones(size(f))];
    beq = load-p_min_total;
    ub = [step1*ones(1,segment_division) step2*ones(1,segment_division) step3*ones(1,segment_division)];
    lb = zeros((size(f)));
    [sol, MinimumCost] = linprog(f, [], [], Aeq, beq, lb, ub);
    
    P1=0; % get final value of three generations 
    P2=0;
    P3=0;
    for i=1:segment_division
        P1=P1+sol(i,1);
        P2=P2+sol(i+segment_division,1);
        P3=P3+sol(i+2*segment_division,1);
    end
    
    P1=P1+p1_min;
    P2=P2+p2_min;
    P3=P3+p3_min;
    
    %Calculate the minimum cost
    syms p1 p2 p3
    F1 = 561 + 7.92 * p1 + 0.001562 * p1^2; %cost funciton1
    F2 = 310 + 7.85 * p2 + 0.00194 *  p2^2; %cost function2
    F3 = 93.6 + 9.564* p3 + 0.005784 *  p3^2;  % cost function3
    
    Minimumcost = double(subs(F1, p1, P1)) + double(subs(F2, p2, P2)) + double(subs(F3, p3, P3));
    %display the final results
     MaxGen=p1_max +  p2_max + p3_max;
    disp(['Max Gen: ' num2str(MaxGen) ' MW']);
      disp(['Min Gen: ' num2str(p_min_total) ' MW']);
    disp(['Optimal P1: ' num2str(P1) ' MW']);
    disp(['Optimal P2: ' num2str(P2) ' MW']);
    disp(['Optimal P3: ' num2str(P3) ' MW']);
      disp(['F1: ' num2str(double(subs(F1, p1, P1))) ' per hour']);
    disp(['F2: ' num2str(double(subs(F2, p2, P2))) ' per hour']);
    disp(['F3: ' num2str(double(subs(F3, p3, P3))) ' per hour']);
    disp(['Minimum Cost: $' num2str(Minimumcost) ' per hour']);
end
    
