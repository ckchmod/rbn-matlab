%This function is a simple template for the data runs we'll take, and is an
%example of how to use the "insert_mutant" function

%% Set up a couple different grids
%   Note: these are random grids

%Original grid
p=0.5;
gridOrig = boolCellGrid('orthogonal',4^2,18,2,p,2,[],[],[]);

%Grid with mutant cells; note that the truth tables and varF will be
%different here, which probably isn't what we want
p=0.7;
gridMut = boolCellGrid('orthogonal',4^2,18,2,p,2,[],[],[]);


%==========================================================================




%% Create grids with different numbers of mutants

%With one mutant
cellPos = 1;
gridOrig_mut1 = copy(gridOrig);
gridOrig_mut1 = gridOrig_mut1.insert_mutants(gridMut,cellPos); %Pass the full grid of mutant cells

%Two mutants
cellPos = [3, 14];
gridOrig_mut2 = copy(gridOrig);
gridOrig_mut2 = gridOrig_mut2.insert_mutants(gridMut,cellPos); %Pass the full grid of mutant cells


%==========================================================================




%% Simulate all the grids

gridOrig.update_all(50);

gridOrig_mut1.update_all(50);

gridOrig_mut2.update_all(50);


%==========================================================================



%% Plot all of them

gridOrig.plot_cells;

gridOrig_mut1.plot_cells;

gridOrig_mut2.plot_cells;

%==========================================================================



