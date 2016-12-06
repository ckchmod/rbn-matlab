%% Mutant Cell Count Trials

%define baseline parameters
iterations = 50;
numcells = 8^2;
numgenes = 8;
k = 1;
p = 0.5;
bandwidth = 5;

%generate homogeneous grid as baseline
Orig = boolCellGrid('orthogonal',numcells,numgenes,k,p,bandwidth,[],[],[]);
Origdummy = copy(Orig);
Origdummy.update_all(iterations); % update original grid

Origdummy.plot_cells
%generate grid of 'mutant cells'
%make some alteration


pchangedistances = zeros(numcells,4);
%create copies of baseline grid and insert appropriate number of mutant cells
% ps = 0.5:0.1:0.8
for m = 1:4
    k = m;
    fprintf('Current k value: %d\n',k);
    gridMut = boolCellGrid('orthogonal',numcells,numgenes,k,p,bandwidth,[],[],[]);
    mutant_cells = cell(numcells,1);
    distances = zeros(numcells,1);
    
    for j = 1:numcells
        cellPos = randi([1 numcells],1,j);
        mutant = copy(Orig);
        mutant_cells{j} = mutant.insert_mutants(gridMut,cellPos); %Pass the full grid of mutant cells
        mutant = mutant_cells{j};
        mutant.update_all(iterations); %update all the copies
        distances(j) = distance(Origdummy,mutant); %calculate hamming distances
    end
    pchangedistances(:,m) = distances;
end


plot(pchangedistances, '*')
title('Ordered GRN with Orthogonal Communication')
xlabel('Number of Mutant Cells')
ylabel('Hamming Distance')
legend('k=1', 'k=2','k=3','k=4')