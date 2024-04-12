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

public class DormantPercentage extends CalculationNode implements Loggable, Function {
    final public Input<SeedbankTree> seedbankTreeInput = new Input<>("sbTree", "tree to report dormant percentage.", Validate.REQUIRED);
    final public Input<Boolean> logLengthsInput = new Input<>("logLengths", 
    		"whether to log lengths of dormant and active sections (default false)", false);
  
    
    // dormant branch length proportion of total tree
    double percentage;
    
    // lengths[0] -- dormant; lengths[1] -- active
    double[] lengths;
    
    boolean logLengths;

	@Override
	public void initAndValidate() {
	    double percentage = 0.0;
		lengths = new double[]{0.0, 0.0};
		logLengths = logLengthsInput.get();
	}
	
	public void update() {
		final SeedbankTree sbTree = seedbankTreeInput.get();
		lengths[0] = 0.0;
		lengths[1] = 0.0;

        for (Node node : sbTree.getNodesAsArray()) {
            if (node.isRoot()) {
                continue;
            }

            SeedbankNode sbNode = (SeedbankNode)node;
            int thisType = sbNode.getNodeType();
            double lastTime = sbNode.getHeight();
            for (int i = 0; i < sbNode.getChangeCount(); i++) {
                double nextTime = sbNode.getChangeTime(i);
                lengths[thisType] += (nextTime - lastTime);
                lastTime = nextTime;
                thisType = 1 - thisType;
            }

            lengths[thisType] += sbNode.getParent().getHeight() - lastTime;
        }
        
        percentage = lengths[0] / (lengths[0] + lengths[1]);
	}
	
	@Override
	public void init(PrintStream out) {
        final Tree tree = seedbankTreeInput.get();
        
        if (logLengths) {
        	out.print(tree.getID() + ".dormantLength\t");
            out.print(tree.getID() + ".activeLength\t");
        }
        
        out.print(tree.getID() + ".dormantPercentage\t");
	}

	@Override
	public void log(long sample, PrintStream out) {
        update();
        
        if (logLengths) {
        	out.print(lengths[0] + "\t");
            out.print(lengths[1] + "\t");
        }
        
        out.print(percentage + "\t");
	}

	@Override
	public void close(PrintStream out) {

	}

	// function implementation
	
	@Override
	public int getDimension() {
		return 1;
	}

	@Override
	public double getArrayValue(int dim) {
		update();
		
		return percentage;
	}

}
