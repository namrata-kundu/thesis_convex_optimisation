function  visualize(f,f_right_bound, rho, rho_right_bound)
    piecewise_f = build_piecewise_function(f, f_right_bound);
    fplot(piecewise_f)
    grid on;
    hold on;
    piecewise_rho = build_piecewise_function(rho, rho_right_bound);

    fplot(piecewise_rho)
end

