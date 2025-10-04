%% Objective Function
function cost = cc_update_cost(norminal_cc, exp_cc, n_select_cc)

    cost = mean((norminal_cc(2:1+n_select_cc) - exp_cc(2:1+n_select_cc)).^2, "all");
    
end
