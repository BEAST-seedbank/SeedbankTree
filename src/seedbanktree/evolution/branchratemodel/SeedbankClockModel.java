package seedbanktree.evolution.branchratemodel;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import seedbanktree.evolution.tree.SeedbankTree;

public class SeedbankClockModel extends BranchRateModel.Base {

	final public Input<Function> alphaInput = new Input<>("alpha", "Dormant rate multiplier");
	final public Input<Function> lambdasInput = new Input<>("lambdas", "Branch dormant fraction");
	final public Input<IntegerParameter> indicatorsInput = new Input<>("indicators", "Dormancy indicators");
	final public Input<SeedbankTree> treeInput =
            new Input<>("tree", "the tree this relaxed clock is associated with.", Validate.REQUIRED);
	
	SeedbankTree tree;
    Function meanRate;
    
	@Override
	public void initAndValidate() {
		tree = treeInput.get();
		
		meanRate = meanRateInput.get();
        if (meanRate == null) {
            meanRate = new RealParameter("1.0");
        }
	}
	
	@Override
	public double getRateForBranch(Node node) {
		// this must be synchronized to avoid being called simultaneously by
        // two different likelihood threads
//    	synchronized (this) {
//    		if (recompute) {
//                recalculateScaleFactor();
//                recompute = false;
//			}
//        }
		
		// Handle root
		if (node.isRoot()) {
			return 1.0;
		}
		
		// Indicators off
		if (indicatorsInput.get().getValue(node.getNr()) == 0) {
			return 1.0;
		}
		
		double lambda = lambdasInput.get().getArrayValue(node.getNr());
		double r = 1 - (1 - alphaInput.get().getArrayValue()) * lambda;
		
        return r;
	}

}
