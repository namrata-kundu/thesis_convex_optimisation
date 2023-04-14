function  [g, g_pieces, objective] = get_nearest_convex_function_with_given_number_of_pieces(f, pieces, num_of_pieces)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%     yalmip('clear')
 
    boundary_limits = 10;
    g = [];
    g_pieces = [];
    repeating_g = [];
    
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

    bounds = [pieces(1,:)];
    bounds(1) = left_boundary;
    bounds(end) = right_boundary;

    total_num_of_pieces_in_f = size(pieces,2);
    each_piece_in_g_equivalent_to_num_of_f_pieces = floor(total_num_of_pieces_in_f/num_of_pieces);

    last_i=0;
    g_pieces = [g_pieces, bounds(1)];
    for i=1+each_piece_in_g_equivalent_to_num_of_f_pieces:each_piece_in_g_equivalent_to_num_of_f_pieces:size(bounds,2)-1
        g_pieces = [g_pieces, bounds(i)];
        a=sdpvar(1,1);
        b=sdpvar(1,1);
        c=sdpvar(1,1);
        col = [a;b;c];
        g = [g, col];
        for j = i:i+each_piece_in_g_equivalent_to_num_of_f_pieces-1
            repeating_g = [repeating_g, col];
        end
        last_i = i;
    end
    
    g_pieces = [g_pieces, bounds(end)]; 
    a=sdpvar(1,1);
    b=sdpvar(1,1);
    c=sdpvar(1,1);
    col = [a;b;c];
     g = [g, col];
    for k = last_i:size(bounds,2)
        repeating_g = [repeating_g, col];
    end

    objective = 0;
    syms a_var b_var c_var x 
    for j=1:each_piece_in_g_equivalent_to_num_of_f_pieces:size(bounds,2)-1 
        for i=j:j+each_piece_in_g_equivalent_to_num_of_f_pieces-1  
            if i>size(bounds,2)-1
                break; %added check so that if i doesn't go beyond range of f
            end
            af=f(1,i);
            bf=f(2,i);
            cf=f(3,i);
            func = ((af*x*x + bf*x + cf) - (a_var*x*x + b_var*x + c_var))^2;
            lower_bound = bounds(i);
            upper_bound = bounds(i+1);
            symbolic_integral = int(func, [lower_bound upper_bound]);
            str_symbolic_integral = char(symbolic_integral);
            str_symbolic_integral = strrep(str_symbolic_integral, 'a_var', strcat('repeating_g(1,',num2str(i),')'));
            str_symbolic_integral = strrep(str_symbolic_integral, 'b_var', strcat('repeating_g(2,',num2str(i),')'));
            str_symbolic_integral = strrep(str_symbolic_integral, 'c_var', strcat('repeating_g(3,',num2str(i),')'));
            integral = eval(str_symbolic_integral);
            objective = objective + integral;
        end
    end

    %Build constraints
    Constraints = [];
    
    %a>=0
    for i=1:size(g,2)
        Constraints = [Constraints, g(1,i)>=0];
    end

    %rho_i(x_(i+1))=rho_(i+1)(x_(i+1))
    for i=1:size(g_pieces,2)-2 %check 1 or 2
        x_val = g_pieces(i+1);
        ai_val = g(1,i);
        bi_val = g(2,i);
        ci_val = g(3,i);
        aiplus1_val = g(1,i+1);
        biplus1_val = g(2,i+1);
        ciplus1_val = g(3,i+1);
        
        Constraints = [Constraints, ai_val*x_val*x_val + bi_val*x_val + ci_val == aiplus1_val*x_val*x_val + biplus1_val*x_val + ciplus1_val ];
    end
 
    % rho_i'(x_(i+1))<=rho_(i+1)'(x_(i+1))
    for i=1:size(g_pieces,2)-2 %check 1 or 2
        x_val = g_pieces(i+1);
        ai_val = g(1,i);
        bi_val = g(2,i);
        aiplus1_val = g(1,i+1);
        biplus1_val = g(2,i+1);
        
        Constraints = [Constraints, 2*ai_val*x_val + bi_val <= 2*aiplus1_val*x_val + biplus1_val ];
    end

    %last pieces should have function same as boundary if boundaries go to
    %infinity
    if isinf(pieces(1))
        constraint1 = [g(1,1)==f(1,1), g(2,1)==f(2,1), g(3,1)==f(3,1)];
        Constraints = [Constraints constraint1];
    end
    if isinf(pieces(end))
        constraint2 = [g(1,end)==f(1,end), g(2,end)==f(2,end), g(3,end)==f(3,end)];
        Constraints = [Constraints constraint2];
    end

    options = sdpsettings('solver','baron', 'verbose',1);
    sol = optimize(Constraints,objective, options)
    disp(value(objective));
    objective = value(objective);

    g = convert_to_values_rho(g);

%     visualize(f,pieces,g,g_pieces);

end




function val = convert_to_values_rho(rho)
    val = [];
    for i=1:size(rho,1)
        val_row = [];
        for j=1:size(rho,2)
            val_row = [val_row value(rho(i,j))];
        end
        val = [val; val_row];
    end        
end

