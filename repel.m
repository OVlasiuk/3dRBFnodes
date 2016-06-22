function cnf = repel(cnf)

k_value = 15;           % number of nearest neighbors used in the knnsearch
repel_steps = 1;

[IDX] = knnsearch(cnf, cnf, 'k', k_value+1);
forces = zeros(size(cnf,1), size(cnf,2));        


% force_i = zeros(k_value, 3);              


for iter=1:repel_steps
    for i=1:size(cnf,1)
        force_i = normr(ones(k_value,1) * cnf(i,:) - cnf(IDX(i,2:end),:));
        factor = 1/(iter*(20 + density(cnf(i,:))));
        forces(i, :) = factor * sum(force_i, 1);
    end
    cnf = cnf + forces;
    forces(4000:4010,:)
end
