package seedbanktree.util;

import java.io.PrintStream;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;
import seedbanktree.evolution.tree.SeedbankNode;
import seedbanktree.evolution.tree.SeedbankTree;

public class NodeTypeCountsLogger extends CalculationNode implements Loggable, Function {
    final public Input<SeedbankTree> seedbankTreeInput = new Input<>("sbTree", "tree to report type lengths for.", Validate.REQUIRED);
    
    // counts[0] = dormant; counts[1] = active
    int[] counts;

	@Override
	public void initAndValidate() {
		counts = new int[]{0, 0};
	}
	
	public void update() {
		final SeedbankTree sbTree = seedbankTreeInput.get();
		counts[0] = 0;
		counts[1] = 0;
		
        for (Node node : sbTree.getNodesAsArray()) {
        	SeedbankNode sbNode = (SeedbankNode)node;
        	counts[sbNode.getNodeType()] += 1;
        			
            if (node.isRoot()) {
                continue;
            }
            
            if (sbNode.getNodeType() == 0) {
//            	assert(sbNode.getChangeCount() % 2 == 1);
            	counts[0] += sbNode.getChangeCount() / 2;
    			counts[1] += 1 + sbNode.getChangeCount() / 2;
            } else {
//            	assert(sbNode.getChangeCount() % 2 == 0);
    			counts[0] += sbNode.getChangeCount() / 2;
    			counts[1] += sbNode.getChangeCount() / 2;
            }
        }
	}
	
	@Override
	public void init(PrintStream out) {
        final Tree tree = seedbankTreeInput.get();
        out.print(tree.getID() + ".dormantNodeCount\t");
        out.print(tree.getID() + ".activeNodeCount\t");
	}

	@Override
	public void log(long sample, PrintStream out) {
        update();
        
        out.print(counts[0] + "\t");
        out.print(counts[1] + "\t");
	}

	@Override
	public void close(PrintStream out) {

	}

	// function implementation
	
	@Override
	public int getDimension() {
		return 2;
	}

	@Override
	public double getArrayValue(int dim) {
		update();
		
		if (dim == 0) 
			return counts[0];
		else if (dim == 1)
			return counts[1];
		else
			throw new IllegalArgumentException("provided dim is " + dim + ", it should be 0 or 1");
	}

}
