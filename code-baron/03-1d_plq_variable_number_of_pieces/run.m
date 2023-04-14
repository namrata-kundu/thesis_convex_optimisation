yalmip('clear')
% W function
pieces = [-inf(), -5, 0, 5, inf()];
f = [0, 0, 0, 0;
    -1, 1, -1, 1;
    -5, 5, 5, -5];
epsilon=0.000005;
[rho,new_pieces,  objective]  = nearest_convex_function_variable_pieces(f,pieces);
[g,g_pieces, num_of_pieces] = decrease_pieces_of_convex_function(rho,new_pieces, epsilon);

% disp(value(objective))
disp(size(new_pieces,2)-1)
disp( num_of_pieces )


g_pieces = convert_to_values(g_pieces);
rho_vals = convert_to_values_rho(g);

visualize(rho,new_pieces,g,g_pieces);


function val = convert_to_values(new_pieces)
    val = [];
    for i=1:size(new_pieces,2)
        val = [val value(new_pieces(i))];
    end        
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
