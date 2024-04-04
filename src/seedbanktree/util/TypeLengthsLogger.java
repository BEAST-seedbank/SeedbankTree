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

public class TypeLengthsLogger extends CalculationNode implements Loggable, Function {
    final public Input<SeedbankTree> seedbankTreeInput = new Input<>("sbTree", "tree to report type lengths for.", Validate.REQUIRED);
    
    // lengths[0] -- dormant; lengths[1] -- active
    double[] lengths;

	@Override
	public void initAndValidate() {
		lengths = new double[]{0.0, 0.0};
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
	}
	
	@Override
	public void init(PrintStream out) {
        final Tree tree = seedbankTreeInput.get();
        out.print(tree.getID() + ".dormantLength\t");
        out.print(tree.getID() + ".activeLength\t");
	}

	@Override
	public void log(long sample, PrintStream out) {
        update();
        
        out.print(lengths[0] + "\t");
        out.print(lengths[1] + "\t");
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
			return lengths[0];
		else if (dim == 1)
			return lengths[1];
		else
			throw new IllegalArgumentException("provided dim is " + dim + ", it should be 0 or 1");
	}

}
