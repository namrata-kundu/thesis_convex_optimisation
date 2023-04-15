function [rho,new_pieces, num_of_pieces] = decrease_pieces_of_convex_function(f,pieces, epsilon)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    r = size(pieces, 2);
%     [rho,new_pieces, num_of_pieces] = binary_search(f, pieces, 1, r, epsilon);
[rho,new_pieces, num_of_pieces] = linear_search(f, pieces, 1, r, epsilon);
    
end

function [rho,new_pieces, mid] = binary_search(f, pieces, l, r, epsilon)

    if r>=l
        mid = l + floor((r-l)/10);
        
        [rho,new_pieces, objective] = get_nearest_convex_function_with_given_number_of_pieces(f, pieces, mid);
        
        if objective <= epsilon
            %rho, new_pieces is our answer
%             [rho,new_pieces] = binary_search(f, pieces, 1, mid-1, epsilon);
        else
             [rho,new_pieces] = binary_search(f, pieces, mid+1, r, epsilon);
        end
    else
        rho = f;
        new_pieces = pieces;
    end
end




function [rho,new_pieces, mid] = linear_search(f, pieces, l, r, epsilon)
    if r>=l
        for i=l:r
            mid=i;
            [rho,new_pieces, objective] = get_nearest_convex_function_with_given_number_of_pieces(f, pieces, mid);
            
            if objective <= epsilon
                break;
            end
        end
    else
        rho = f;
        new_pieces = pieces;
    end
end
    