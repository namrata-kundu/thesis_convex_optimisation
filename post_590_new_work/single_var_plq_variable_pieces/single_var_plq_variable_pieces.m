

yalmip('clear')
 
%%%%%%%start here

boundary_limits = 10;
upper_boundary = inf();

%W function
% pieces = [-inf(), -5, 0, 5, inf()];
% f = [0, 0, 0, 0;
%     -1, 1, -1, 1;
%     -5, 5, 5, -5];

%absolute value function
% pieces = [-inf(), 0, inf()];
% f = [
%     0, 0;
%     -1, 1;
%     0, 0];
% 
% %convex function
pieces = [-inf(), 1, inf()];
f = [
    (1/2), 0;
    0, 2;
    1, -(1/2)];
 

% %linear function
pieces = [-inf(), 0, inf()];
f = [0, 0;
    1, 1;
    0, 0];
% 
% 
% %linear function -- doesn't work right now
% f = [-inf();
%     0;
%     1;
%     0];


% % %the very first function I started with
% pieces = [-inf(), 1, 2.5, 6, inf()];
% f = [
%     (1/2), 0, 0, 0;
%     0, 2, -1, 1;
%     1, -(1/2), 7, -5];
% 

% %the very first function I started with
pieces = [-inf(), 1, 2.5, 6];
f = [
    (1/2), 0, 0;
    0, 2, -1;
    1, -(1/2), 7];


[row_size,col_size] = size(f);

%Boundary pieces convexity check
if f(1,1)<0 || f(1,col_size)<0 %a1<0, an<0
    disp("Infeasible problem: Boundary pieces are not convex")
    return;
end

if size(f,2)==1
    %check what needs to be done
end

%check if left tangent > right tangent
if isinf(pieces(1)) && isinf(pieces(end))
    if 2*f(1,1)*pieces(2) + f(2,1) > 2*f(1,col_size)*pieces(end) + f(2,col_size)
        disp("Infeasible problem: left boundary piece tangent > right boundary piece tangent")
        return;
    end
end

number_of_function_pieces = col_size;
new_number_of_function_pieces = number_of_function_pieces*10; %each piece will be divided into 10 pieces

rho=[]; %resultant function matrix
new_f = []; %(f divided into multiple pieces)

% boundary_limits = f(1,col_size) - f(1,2); %f(1,col_size) will always be > f(1,2)
if pieces(1) == -inf()
    left_boundary = pieces(2) - boundary_limits; 
else
    left_boundary = pieces(1);
end

if pieces(end) == inf()
    right_boundary = pieces(end-1) + boundary_limits;
else
    right_boundary = pieces(end);
end
    

% rho = [rho, linspace(left_boundary,f(1,2),round(boundary_limits))];
bounds = [pieces(1,:)];
bounds(1) = left_boundary;
bounds(end) = right_boundary;

new_pieces = [];
granularity_of_pieces = 1; %(default to 1 - for having piece every 1 point)
% construct new_f (f divided into multiple pieces)
for i=1:size(bounds,2)-1
    left_bound = bounds(i);
    right_bound = bounds(i+1)-0.005; %Added 0.005 to not have repetative boundary points
    num_of_pieces = round(right_bound - left_bound)*granularity_of_pieces;
    first_bound_row = linspace(left_bound,right_bound,(num_of_pieces));
    a_row = [f(1,i)*ones(1,(num_of_pieces))]; 
    b_row = [f(2,i)*ones(1,(num_of_pieces))]; 
    c_row = [f(3,i)*ones(1,(num_of_pieces))]; 
    new_pieces = [new_pieces first_bound_row];
%     values = [first_bound_row; a_row; b_row; c_row];
    values = [a_row; b_row; c_row];
    new_f = [new_f values];
end

new_pieces = [new_pieces right_boundary];


% Allocate sdpvar variables
a = cell(size(new_f,2),1);
b = cell(size(new_f,2),1);
c = cell(size(new_f,2),1);
for i=1:size(new_f,2)
%     first_bound_row = new_f(1,i);
    a{i} = sdpvar(1,1); 
    b{i} = sdpvar(1,1); 
    c{i} = sdpvar(1,1); 
%     values = [first_bound_row; a{i}; b{i}; c{i}];
    values = [a{i}; b{i}; c{i}];
    rho = [rho values];
end
%%Remember - last column value is extra under the right bound (not needed)


%Find integrals - area using symbolic integrals
objective = 0;
syms a_var b_var c_var x
for i=1:size(new_f,2) %total_number_of_pieces = size(new_f,2);
    af=new_f(1,i);
    bf=new_f(2,i);
    cf=new_f(3,i);
    func = ((af*x*x + bf*x + cf) - (a_var*x*x + b_var*x + c_var))^2;
    lower_bound = new_pieces(i);
    upper_bound = new_pieces(i+1);
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
for i=1:size(new_f,2)
    Constraints = [Constraints, rho(1,i)>=0];
end

 
%rho_i(x_(i+1))=rho_(i+1)(x_(i+1))
% for i=1:size(new_f,2)-2
for i=1:size(new_f,2)-1
    x_val = new_pieces(i+1);
    ai_val = rho(1,i);
    bi_val = rho(2,i);
    ci_val = rho(3,i);
    aiplus1_val = rho(1,i+1);
    biplus1_val = rho(2,i+1);
    ciplus1_val = rho(3,i+1);
    
    Constraints = [Constraints, ai_val*x_val*x_val + bi_val*x_val + ci_val == aiplus1_val*x_val*x_val + biplus1_val*x_val + ciplus1_val ];
end



 
% rho_i'(x_(i+1))<=rho_(i+1)'(x_(i+1))
for i=1:size(new_f,2)-1
    x_val = new_pieces(i+1);
    ai_val = rho(1,i);
    bi_val = rho(2,i);
    aiplus1_val = rho(1,i+1);
    biplus1_val = rho(2,i+1);
    
    Constraints = [Constraints, 2*ai_val*x_val + bi_val <= 2*aiplus1_val*x_val + biplus1_val ];
end


%last pieces should have function same as boundary
if isinf(pieces(1))
    constraint1 = [rho(1,1)==f(1,1), rho(2,1)==f(2,1), rho(3,1)==f(3,1)];
    Constraints = [Constraints constraint1];
end
if isinf(pieces(end))
    constraint2 = [rho(1,end)==f(1,end), rho(2,end)==f(2,end), rho(3,end)==f(3,end)];
    Constraints = [Constraints constraint2];
end 





%Feed into solver
options = sdpsettings('solver','', 'verbose', 1);
sol = optimize(Constraints,objective, options)
disp(value(objective))
% visualize(f,bounds(size(bounds,2)),rho,rho(1,size(rho,2)));
visualize(f,pieces,rho,new_pieces);




