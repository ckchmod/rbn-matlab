Readme for the boolCellGrid.m grid class and the boolCell.m class

12/1/2016

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


Copy and paste the MATLAB code above to run some examples!

A couple of the choices I've made that may or may not be reasonable:
	The Ttable is randomly filled, which means genes are turned 'on' rather more than in the original network
	

TO DO:
	Implement other topologies, which can be easily added to the 'setNeighbors' function in the file 'boolCellGrid'
		This function is explained in its header comment, but I'd love to talk about it more

	Prettier visualization; easiest way is to just add a 'plot_cells2' function to the 'boolCellGrid' class

	Actually take data and make sense of it!
		Two nice features: 
			Each of the cells themselves have properties that can be changed, which can be done in the middle of a simulation
			When you run "a.update_all(x)" you update the grid cells for 'x' timesteps and it keeps track of the current time.
			Thus we can just 'update_all' for a certain number of steps, change a cell's properties, and continue to simulate
			