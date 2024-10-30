package seedbanktree.evolution.branchratemodel;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import seedbanktree.evolution.tree.SeedbankTree;

@Description("A branch rate clock model for seedbank trees. It calculates the local branch rate multiplier "
		+ "based on the amount of dormancy on a branch (parameterized by a spike and slab prior) and the "
		+ "dormant mutation rate multiplier.")
public class SeedbankClockModelSNS extends BranchRateModel.Base {

	final public Input<Function> alphaInput = 
			new Input<>("alpha", "Dormant mutation rate multiplier", Validate.REQUIRED);
	
	final public Input<RealParameter> lambdasInput = 
			new Input<>("lambdas", "Branch dormant fraction", Validate.REQUIRED);
	
	final public Input<IntegerParameter> etasInput = 
			new Input<>("etas", "Spike and slab mixture indicator", Validate.REQUIRED);
	
	final public Input<SeedbankTree> sbTreeInput = 
			new Input<>("sbTree", "the seedbank tree this relaxed clock is associated with.", Validate.REQUIRED);
	
	SeedbankTree sbTree;
    Function meanRate;
    
	@Override
	public void initAndValidate() {
		sbTree = sbTreeInput.get();
		
		meanRate = meanRateInput.get();
        if (meanRate == null) {
            meanRate = new RealParameter("1.0");
        }
	}
	
	@Override
	public double getRateForBranch(Node node) {
		// Handle root
		if (node.isRoot()) {
			return 1.0;
		}
		
		// If eta is 0, the local branch rate is 1
		if (etasInput.get().getValue(node.getNr()) == 0) {
			return 1.0;
		}
		
		// If eta is 1, calculate the local branch rate
		double lambda = lambdasInput.get().getArrayValue(node.getNr());
		double r = 1 - (1 - alphaInput.get().getArrayValue()) * lambda;
		
        return r;
	}

}
