Readme for the boolCellGrid.m grid class and the boolCell.m class

LAST UPDATE: 12/7/2016

Here are a couple examples to get everybody started, and they're at the top of the boolCellGrid class:
    
    EXAMPLE1 - 
    % Create a RBN
    
    a=boolCellGrid('symmetric',4,18,2,.5,1, [], [], []); 
    a.update_all(50); 
    a.plot_cells;
     

    EXAMPLE2 -
    % Store arrays (Initial States, Truth Table, Wiring Nodes) and run the specified parameters
    
    i = a.initStates;
    t = a.initTtable;
    v = a.initvarF;
    b = boolCellGrid('symmetric',4,18,2,.5,1, i, t, v);

    EXAMPLE3 -
    % Saves movie at FPS = 1/dt.

    b.plotCells(true, .5);

    EXAMPLE4 - 
    % Create a RBN with perturbation (small probability of gene flipping at each iteraiton) and find the Steady State Distribution of each cell.
    % Default perturbation (no parameter specfied) is 0.

    pert = .2;
    a=boolCellGrid('symmetric',4,18,2,.5,1, [], [], [], pert); 
    a.update_all(50);
    ssDist(a) 

Copy and paste the MATLAB code above to run some examples!

TO DO:

	Lypanuov Exponent

	Different Communication Function

	Diversity

Two nice features: 

	Each of the cells themselves have properties that can be changed, which can be done in the middle of a simulation

	When you run "a.update_all(x)" you update the grid cells for 'x' timesteps and it keeps track of the current time.

	Thus we can just 'update_all' for a certain number of steps, change a cell's properties, and continue to simulate


