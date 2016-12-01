Readme for the boolCellGrid.m grid class and the boolCell.m class

11/29/2016

Here are a couple examples to get everybody started, and they're at the top of the boolCellGrid class:

    
      EXAMPLE2 - 
	a=boolCellGrid('line',4,18,2,.5,1, [], [], []); a.update_all(50); a.plot_cells(1.0);
              This random network seems to produce an oscillation.
		The optional argument in the method 'plot_cells' will pause the graph

      EXAMPLE3 -
	a=boolCellGrid('symmetric',5^2,18,2,3,false,42); a.update_all(20); a.plot_cells;
	      This is a matrix of cells (5x5) that is visualized as we mentioned before class, 
	    	with a unique color for each cell state. It's not as pretty as we might want it though,
		and the states can blend into each other pretty easily


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
			