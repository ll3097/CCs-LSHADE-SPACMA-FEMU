clear all
%% Main program for LSHADE-SPACMA-CC FEMU
% This script implements the customized cepstrum-based LSHADE-SPACMA hybrid FEMU framework 
% for structural parameter identification problems in SHM.

%% Load experimental data (FEMU objective)
% Example: baseline acceleration response from LANL three-story frame dataset.
load('exp_resp_baseline.mat')
exp_cc = extract_cc(exp_resp);   % Extract cepstral coefficients from experimental response

%% Structural parameter configuration
n_select_cc = 100;       % Number of cepstral coefficients to use in the cost function
DOF = 8;                 % Degrees of freedom of the structure
dim = DOF*2;             % Optimization dimension (stiffness + damping for each DOF)

% Parameter bound settings for random initialization and optimization process 
% e.g., ±20-40% for stiffness
k0 = 25000;              
xi0 = 0.01;
k = repmat(k0, 1, DOF);   % stiffness values (N/m) for 8 DOFs
xi = repmat(xi0, 1, DOF); % damping ratios for 8 DOFs
lb = [k * 0.7, xi * 0.99];  
ub = [k * 1.3, xi * 1.01];  

% % Uncomment to generate additional customized/simulated experimental data instead:
% exp_para = [k, xi];
% exp_resp = response_simulation8d(exp_para, DOF);
% exp_cc = extract_cc(exp_resp);

%% LSHADE-SPACMA parameter settings
pop_size = 50;       % Population size
p_best_rate = 0.2;   % Top p-best rate for mutation
arc_rate = 1.4;      % Archive size rate for external population
memory_size = 10;    % Historical memory size for control parameters
max_pop_size = pop_size;
min_pop_size = 15.0;
L_Rate = 0.80;       % Learning rate for memory adaptation

%% Hybridization settings (DE + CMA-ES)
First_class_percentage = 0.5; % Initial probability for selecting DE operators

%% Initialize population
% Randomly initialize candidate solutions within parameter bounds
popold = repmat(lb, pop_size, 1) + rand(pop_size, dim) .* (repmat(ub - lb, pop_size, 1));
pop = popold;

% Evaluate initial fitness values based on cepstral coefficient matching
fitness = zeros(pop_size, 1);
for it = 1:pop_size
    candidate = pop(it, :);
    sim_resp = response_simulation8d(candidate, DOF);
    sim_cc = extract_cc(sim_resp);
    fitness(it) = cc_update_cost(sim_cc, exp_cc, n_select_cc);
end

% Initialize best-so-far (bsf) solution
bsf_fit_var = inf;         % Best objective value found
bsf_solution = zeros(1, dim);
bsf_index = 0;

for i = 1:pop_size
    if (fitness(i) < bsf_fit_var && isreal(pop(i, :)) && all(~isnan(pop(i, :))))
        bsf_fit_var = fitness(i);
        bsf_solution = pop(i, :);
        bsf_index = i;
    end
end

%% Initialize LSHADE-SPACMA memory and archive
memory_sf = ones(memory_size, 1) .* 1.0;   % Memory of scaling factors
memory_cr = ones(memory_size, 1) .* 0.9;    % Memory of crossover rates
memory_pos = 1;

archive.NP = arc_rate * pop_size;   % Maximum archive size
archive.pop = zeros(0, dim);        % Archived solutions
archive.funvalues = zeros(0, 1);    % Objective values of archive

memory_1st_class_percentage = First_class_percentage .* ones(memory_size, 1);

%% Initialize CMA-ES parameters
sigma = 0.5;           % Step size (standard deviation)
xmean = rand(dim, 1);  % Initial mean point
mu = floor(pop_size/2);

% Recombination weights
weights = log(mu+0.5) - log(1:mu)';
weights = weights / sum(weights);
mueff = sum(weights)^2 / sum(weights.^2); % Effective sample size

% Strategy parameter settings for CMA-ES adaptation
cc = (4 + mueff/dim) / (dim + 4 + 2*mueff/dim); 
cs = (mueff+2) / (dim+mueff+5);  
c1 = 2 / ((dim+1.3)^2 + mueff);  
cmu = min(1-c1, 2 * (mueff-2+1/mueff) / ((dim+2)^2+mueff));
damps = 1 + 2*max(0, sqrt((mueff-1)/(dim+1))-1) + cs;

% Internal CMA-ES variables
pc = zeros(dim, 1); ps = zeros(dim, 1);
B = eye(dim); D = ones(dim, 1);
C = B * diag(D.^2) * B';
invsqrtC = B * diag(D.^-1) * B';
eigeneval = 0;
chiN = dim^0.5 * (1 - 1/(4*dim) + 1/(21*dim^2));

%% Iteration settings
max_iter = 70;             % Maximum number of generations
Hybridization_flag = 1;    % Activate DE/CMA-ES hybridization
nfes = 0;                  % Number of function evaluations

%% Main optimization loop
% Each generation updates the population using LSHADE-SPACMA operators 
% (mutation, crossover, selection, memory adaptation, population resizing),
% with optional hybridization using CMA-ES steps.
while nfes < max_iter
    
    %%%% Step 1: Sort population by fitness
    % Keep track of the best-performing solutions
    pop = popold; 
    [temp_fit, sorted_index] = sort(fitness, 'ascend');
    
    %%%% Step 2: Parameter sampling from memory
    % Historical scaling factors (sf) and crossover rates (cr) are sampled 
    % from memory arrays with random indices to maintain diversity.
    mem_rand_index = ceil(memory_size * rand(pop_size, 1));
    mu_sf = memory_sf(mem_rand_index);
    mu_cr = memory_cr(mem_rand_index);
    mem_rand_ratio = rand(pop_size, 1);
    
    % Generate crossover rate (normally distributed around memory values)
    cr = normrnd(mu_cr, 0.1);
    term_pos = find(mu_cr == -1);
    cr(term_pos) = 0;
    cr = min(cr, 1);
    cr = max(cr, 0);
    
    % Generate scaling factor (sf) using two different strategies:
    %   - Early stage: fixed small range for stability
    %   - Later stage: Cauchy distribution for exploration
    if(nfes <= max_iter/2)
        sf=0.45+.1*rand(pop_size, 1);
        pos = find(sf <= 0);
        
        while ~ isempty(pos)
            sf(pos)=0.45+0.1*rand(length(pos), 1);
            pos = find(sf <= 0);
        end
    else
        sf = mu_sf + 0.1 * tan(pi * (rand(pop_size, 1) - 0.5));
        
        pos = find(sf <= 0);
        
        while ~ isempty(pos)
            sf(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
            pos = find(sf <= 0);
        end
    end
    sf = min(sf, 1);
    
    %%%% Step 3: Hybridization class assignment
    % Solutions are probabilistically assigned to DE-based or CMA-ES-based 
    % variation operators. When Hybridization_flag = 0, CMA-ES is disabled.
    Class_Select_Index=(memory_1st_class_percentage(mem_rand_index)>=mem_rand_ratio);
    if(Hybridization_flag==0)
        Class_Select_Index=or(Class_Select_Index,~Class_Select_Index);%All will be in class#1
    end
    
    %%%% Step 4: Mutation and crossover
    % Generate trial vectors (vi) using either:
    %   - Differential Evolution (DE/current-to-pbest/1/bin)
    %   - CMA-ES sampling from a Gaussian distribution
    r0 = 1:pop_size;
    popAll = [pop; archive.pop];
    [r1, r2] = gnR1R2(pop_size, size(popAll, 1), r0);
    
    % Choose top p-best solutions
    pNP = max(round(p_best_rate * pop_size), 2);
    randindex = ceil(rand(1, pop_size) .* pNP);
    pbest = pop(sorted_index(randindex), :);
    
    % --- DE-based offspring ---
    vi=[];
    temp=[];
    if(sum(Class_Select_Index)~=0)
        vi(Class_Select_Index,:)=pop(Class_Select_Index,:)+sf(Class_Select_Index,ones(1, dim)).*... 
        (pbest(Class_Select_Index,:) - pop(Class_Select_Index,:) + pop(r1(Class_Select_Index), :) - popAll(r2(Class_Select_Index), :));
    end
    
    if(sum(~Class_Select_Index)~=0)
        for k=1:sum(~Class_Select_Index)
            temp(:,k) = xmean + sigma * B * (D .* randn(dim,1)); % m + sig * Normal(0,C)
        end
        vi(~Class_Select_Index,:) = temp';
    end
    
    if(~isreal(vi))
        Hybridization_flag=0;
        continue;
    end

    % Ensure within the up and low bound
%     vi = boundConstraint(vi, pop, lu);
    vi_size = size(vi);
    for s = 1:vi_size(1)
        vi(s,:) = max(min(vi(s,:), ub), lb); % This Enforce bound is really important!
    end
    
    % Perform binomial crossover between parent and mutant
    mask = rand(pop_size, dim) > cr(:, ones(1, dim));
    rows = (1 : pop_size)'; 
    cols = floor(rand(pop_size, 1) * dim)+1; 
    jrand = sub2ind([pop_size dim], rows, cols); 
    mask(jrand) = false;
    ui = vi; ui(mask) = pop(mask);
    
    %%%% Step 5: Fitness evaluation of offspring
    children_fitness = zeros(pop_size,1);
    for it = 1:pop_size
        candidate = ui(it,:);
        sim_resp = response_simulation8d(candidate, DOF);
        sim_cc = extract_cc(sim_resp);
        children_fitness(it) = cc_update_cost(sim_cc, exp_cc, n_select_cc);
    end
    
    %%%% Step 6: Update best solution
    % Compare children with best-so-far solution
    for i = 1 : pop_size
        if (children_fitness(i) < bsf_fit_var && isreal(ui(i, :)) && all(~isnan(ui(i, :))))
            bsf_fit_var = children_fitness(i);
            bsf_solution = ui(i, :);
            bsf_index = i;
        end
    end
    
    %%%% Step 7: Selection (survival of the fittest)
    % Compare parent vs. offspring, keep the better one
    dif = abs(fitness - children_fitness);

    Child_is_better_index = (fitness > children_fitness);
    goodCR = cr(Child_is_better_index == 1);
    goodF = sf(Child_is_better_index == 1);
    dif_val = dif(Child_is_better_index == 1);
    dif_val_Class_1 = dif(and(Child_is_better_index,Class_Select_Index) == 1);
    dif_val_Class_2 = dif(and(Child_is_better_index,~Class_Select_Index) == 1);
    
    archive = updateArchive(archive, popold(Child_is_better_index == 1, :), fitness(Child_is_better_index == 1));
    
    [fitness, Child_is_better_index] = min([fitness, children_fitness], [], 2);
    popold = pop;
    popold(Child_is_better_index == 2, :) = ui(Child_is_better_index == 2, :);
    
    num_success_params = numel(goodCR);
    
    %%%% Step 8: Memory update (for sf, cr, and hybridization class probability)
    % Successful parameters are weighted by fitness improvement
    if num_success_params > 0
        sum_dif = sum(dif_val);
        dif_val = dif_val / sum_dif;
        
        % for updating the memory of scaling factor
        memory_sf(memory_pos) = (dif_val' * (goodF .^ 2)) / (dif_val' * goodF);
        
        % for updating the memory of crossover rate
        if max(goodCR) == 0 || memory_cr(memory_pos)  == -1
            memory_cr(memory_pos)  = -1;
        else
            memory_cr(memory_pos) = (dif_val' * (goodCR .^ 2)) / (dif_val' * goodCR);
        end
        
        if (Hybridization_flag==1)% if the Hybridization is activated
            memory_1st_class_percentage(memory_pos) = memory_1st_class_percentage(memory_pos)*L_Rate+ (1-L_Rate)*(sum(dif_val_Class_1) / (sum(dif_val_Class_1) + sum(dif_val_Class_2)));
            memory_1st_class_percentage(memory_pos)=min(memory_1st_class_percentage(memory_pos),0.8);
            memory_1st_class_percentage(memory_pos)=max(memory_1st_class_percentage(memory_pos),0.2);
        end
        
        memory_pos = memory_pos + 1;
        if memory_pos > memory_size;  memory_pos = 1; end
    end
    
    %%%% Step 9: Population size reduction (linear from max_pop_size → min_pop_size)
    % Gradually reduce population size to focus search as iteration progresses
    plan_pop_size = round((((min_pop_size - max_pop_size) / max_iter) * nfes) + max_pop_size);
    
    if pop_size > plan_pop_size
        reduction_ind_num = pop_size - plan_pop_size;
        if pop_size - reduction_ind_num <  min_pop_size; reduction_ind_num = pop_size - min_pop_size;end
        
        pop_size = pop_size - reduction_ind_num;
        for r = 1 : reduction_ind_num
            [valBest, indBest] = sort(fitness, 'ascend');
            worst_ind = indBest(end);
            popold(worst_ind,:) = [];
            pop(worst_ind,:) = [];
            fitness(worst_ind,:) = [];
            Child_is_better_index(worst_ind,:) = [];
        end
        
        archive.NP = round(arc_rate * pop_size);
        
        if size(archive.pop, 1) > archive.NP
            rndpos = randperm(size(archive.pop, 1));
            rndpos = rndpos(1 : archive.NP);
            archive.pop = archive.pop(rndpos, :);
        end

        % update CMA parameters
        mu = pop_size/2;               % number of parents/points for recombination
        weights = log(mu+1/2)-log(1:mu)'; % muXone array for weighted recombination
        mu = floor(mu);
        weights = weights/sum(weights);     % normalize recombination weights array
        mueff=sum(weights)^2/sum(weights.^2); % variance-effectiveness of sum w_i x_i
    end
    
    %%%% Step 10: CMA-ES adaptation (if hybridization active)
    % Update CMA-ES parameters (mean, covariance, step size) 
    % using standard evolution path and recombination rules.
    if(Hybridization_flag==1)
        % Sort by fitness and compute weighted mean into xmean
        [~, popindex] = sort(fitness);  % minimization
        xold = xmean;
        xmean = popold(popindex(1:mu),:)' * weights;  % recombination, new mean value
        
        % Cumulation: Update evolution paths
        ps = (1-cs) * ps ...
            + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean-xold) / sigma;
        hsig = sum(ps.^2)/(1-(1-cs)^(2*nfes/pop_size))/dim < 2 + 4/(dim+1);
        pc = (1-cc) * pc ...
            + hsig * sqrt(cc*(2-cc)*mueff) * (xmean-xold) / sigma;
        
        % Adapt covariance matrix C
        artmp = (1/sigma) * (popold(popindex(1:mu),:)' - repmat(xold,1,mu));  % mu difference vectors
        C = (1-c1-cmu) * C ...                   % regard old matrix
            + c1 * (pc * pc' ...                % plus rank one update
            + (1-hsig) * cc*(2-cc) * C) ... % minor correction if hsig==0
            + cmu * artmp * diag(weights) * artmp'; % plus rank mu update
        
        % Adapt step size sigma
        sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1));
        
        % Update B and D from C
        if nfes - eigeneval > pop_size/(c1+cmu)/dim/10  % to achieve O(problem_size^2)
            eigeneval = nfes;
            C = triu(C) + triu(C,1)'; % enforce symmetry
            if(sum(sum(isnan(C)))>0 || sum(sum(~isfinite(C)))>0 || ~isreal(C))
                Hybridization_flag=0;
                continue;
            end
            [B,D] = eig(C);           % eigen decomposition, B==normalized eigenvectors
            D = sqrt(diag(D));        % D contains standard deviations now
            invsqrtC = B * diag(D.^-1) * B';
        end
        
    end
    
    %%%% Step 11: Record progress
    history(nfes+1) = bsf_fit_var;
    nfes = nfes+1;
    
    fprintf('Iteration %d: Best Cost = %.6f\n', nfes, bsf_fit_var);
    disp(bsf_solution(1:DOF));             % Identified stiffness
    disp(bsf_solution(DOF+1:DOF*2));       % Identified damping ratios
    Convergence_curve(nfes) = bsf_fit_var;
    
end % end main optimization loop


