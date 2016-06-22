function cnf = repel(cnf, r)

k_value = 15;           % number of nearest neighbors used in the knnsearch

[IDX] = knnsearch(transpose(cnf), transpose(cnf), 'k', k_value);
forces = zeros(3, size(cnf,2));

for i=1:length(cnf)
    forces(:, i )
end

