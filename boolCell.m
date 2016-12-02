classdef boolCell < matlab.mixin.Copyable
    %Boolean cell: this object models a cell with a boolean network
    
    properties
        %Set by the caller
        numGenes    %The number of genes in the cell
        k           %The number of connections each node has
        p           %Probability of function value taking on value of 1
        bandwidth   %The number of external connections the cell makes
        Ttable      %The truth table for the network
        varF        %Which intracellular connections a cell makes
        outCells    %Which cells produce intercellular output
        inCells     %Which cells receive intercellular input
        
        %Generated if a grid is set up
        neighbors   %The neighbors of the cell
        cellPos     %Position of the cell
        
        %Generated every step
        states
    end
    
    methods
        
        %Constructor
        function obj = boolCell(numGenes,k,p,bandwidth,Ttable,varF,outCells,inCells)
            %Set some defaults
            if  nargin == 0
                obj.numGenes  = 15;
                obj.k         = 0;
                obj.p         = 0;
                obj.bandwidth = 0;
                obj.Ttable    = [];
                obj.varF      = [];
                obj.outCells  = [];
                obj.inCells   = [];
            else
                obj.numGenes  = numGenes;
                obj.k         = k;
                obj.p         = p;
                obj.bandwidth = bandwidth;
                obj.Ttable    = Ttable;
                obj.varF      = varF;
                obj.outCells  = outCells;
                obj.inCells   = inCells;
            end
        end
        
        function thisCell = setPos(thisCell, cellPos)
            thisCell.cellPos = cellPos;
        end
        
        function setState(thisCell, state, timestep)
            %Used to set the initial state
            thisCell.states(:,timestep) = state;
        end

        function update_genes(thisCell,timestep)
            %-------------- Boolean update for all cells ------------%
            
                for genecol = 1:thisCell.numGenes
                    tempVarf = [];
                    for generow = 1:size(thisCell.varF,1)
                        if (thisCell.varF(generow, genecol) == -1)
                            % do nothing
                        else
                            %This is the list of nodes we will be making
                            %connections to
                            tempVarf = [tempVarf, thisCell.varF(generow,genecol)];
                        end
                    end
                    
                    if (isempty(tempVarf) == 0)
                        wiringnode = 0;
                        for wire = 1:length(tempVarf)
                            % The wiring of input nodes to compute
                            % thisCell's state at t+1
                            wiringnode = wiringnode + 2^(length(tempVarf)-wire) * thisCell.states(tempVarf(wire),timestep-1);
                        end
                        thisCell.states(genecol,timestep) = thisCell.Ttable(wiringnode+1, genecol);
                    end
                end
            
        end
        
    end
end
    
