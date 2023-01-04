

yalmip('clear')
 
%%%%%%%start here

boundary_limits = 10;
upper_boundary = inf();

%W function
f = [-inf(), -5, 0, 5;
    0, 0, 0, 0;
    -1, 1, -1, 1;
    -5, 5, 5, -5];

% %absolute value function
% f = [-inf(), 0;
%     0, 0;
%     -1, 1;
%     0, 0];
% 
% %convex function
% % f = [-inf(), 1;
% %     (1/2), 0;
% %     0, 2;
% %     1, -(1/2)];
% 
% 
% %linear function
% f = [-inf(), 0;
%     0, 0;
%     1, 1;
%     0, 0];
% 
% 
% %linear function -- doesn't work right now
% f = [-inf();
%     0;
%     1;
%     0];


% %the very first function I started with
% f = [-inf(), 1, 2.5, 6;
%     (1/2), 0, 0, 0;
%     0, 2, -1, 1;
%     1, -(1/2), 7, -5];


[row_size,col_size] = size(f);

%Boundary pieces convexity check
if f(2,1)<0 || f(2,col_size)<0 %a1<0, an<0
    disp("Infeasible problem: Boundary pieces are not convex")
    return;
end

if size(f,2)==1
    %check what needs to be done
end

%check if left tangent > right tangent
if 2*f(2,1)*f(1,2) + f(3,1) > 2*f(2,col_size)*f(1,col_size) + f(3,col_size)
    disp("Infeasible problem: left boundary piece tangent > right boundary piece tangent")
    return;
end

number_of_function_pieces = col_size;
new_number_of_function_pieces = number_of_function_pieces*10; %each piece will be divided into 10 pieces

rho=[]; %resultant function matrix
new_f = []; %(f divided into multiple pieces)

% boundary_limits = f(1,col_size) - f(1,2); %f(1,col_size) will always be > f(1,2)
left_boundary = f(1,2) - boundary_limits; 
right_boundary = f(1,col_size) + boundary_limits;

% rho = [rho, linspace(left_boundary,f(1,2),round(boundary_limits))];
bounds = [f(1,:), right_boundary];
bounds(1) = left_boundary;


granularity_of_pieces = 10; %(default to 1 - for having piece every 1 point)
% construct new_f (f divided into multiple pieces)
for i=1:size(bounds,2)-1
    left_bound = bounds(i);
    right_bound = bounds(i+1)-0.005; %Added 0.005 to not have repetative boundary points
    num_of_pieces = round(right_bound - left_bound)*granularity_of_pieces;
    first_bound_row = linspace(left_bound,right_bound,(num_of_pieces));
    a_row = [f(2,i)*ones(1,(num_of_pieces))]; 
    b_row = [f(3,i)*ones(1,(num_of_pieces))]; 
    c_row = [f(4,i)*ones(1,(num_of_pieces))]; 
    values = [first_bound_row; a_row; b_row; c_row];
    new_f = [new_f values];
end


% Allocate sdpvar variables
a = cell(size(new_f,2),1);
b = cell(size(new_f,2),1);
c = cell(size(new_f,2),1);
for i=1:size(new_f,2)
    first_bound_row = new_f(1,i);
    a{i} = sdpvar(1,1); 
    b{i} = sdpvar(1,1); 
    c{i} = sdpvar(1,1); 
    values = [first_bound_row; a{i}; b{i}; c{i}];
    rho = [rho values];
end
%%Remember - last column value is extra under the right bound (not needed)


%Find integrals - area using symbolic integrals
objective = 0;
syms a_var b_var c_var x
for i=1:size(new_f,2)-1 %total_number_of_pieces = size(new_f,2);
    af=new_f(2,i);
    bf=new_f(3,i);
    cf=new_f(4,i);
    func = ((af*x*x + bf*x + cf) - (a_var*x*x + b_var*x + c_var))^2;
    lower_bound = new_f(1,i);
    upper_bound = new_f(1,i+1);
    symbolic_integral = int(func, [lower_bound upper_bound]);
    str_symbolic_integral = char(symbolic_integral);
    str_symbolic_integral = strrep(str_symbolic_integral, 'a_var', strcat('a{',num2str(i),'}'));
    str_symbolic_integral = strrep(str_symbolic_integral, 'b_var', strcat('b{',num2str(i),'}'));
    str_symbolic_integral = strrep(str_symbolic_integral, 'c_var', strcat('c{',num2str(i),'}'));
    integral = eval(str_symbolic_integral);
    objective = objective + integral;
end


%Build constraints
Constraints = [];

%a>=0
for i=1:size(new_f,2)-1
    Constraints = [Constraints, rho(2,i)>=0];
end

 
%rho_i(x_(i+1))=rho_(i+1)(x_(i+1))
for i=1:size(new_f,2)-2
    x_val = rho(1,i+1);
    ai_val = rho(2,i);
    bi_val = rho(3,i);
    ci_val = rho(4,i);
    aiplus1_val = rho(2,i+1);
    biplus1_val = rho(3,i+1);
    ciplus1_val = rho(4,i+1);
    
    Constraints = [Constraints, ai_val*x_val*x_val + bi_val*x_val + ci_val == aiplus1_val*x_val*x_val + biplus1_val*x_val + ciplus1_val ];
end



 
% rho_i'(x_(i+1))<=rho_(i+1)'(x_(i+1))
for i=1:size(new_f,2)-2
    x_val = rho(1,i+1);
    ai_val = rho(2,i);
    bi_val = rho(3,i);
    aiplus1_val = rho(2,i+1);
    biplus1_val = rho(3,i+1);
    
    Constraints = [Constraints, 2*ai_val*x_val + bi_val <= 2*aiplus1_val*x_val + biplus1_val ];
end


%Feed into solver
options = sdpsettings('solver','', 'verbose', 1);
sol = optimize(Constraints,objective, options)
disp(value(objective))
visualize(f,bounds(size(bounds,2)),rho,rho(1,size(rho,2)));




