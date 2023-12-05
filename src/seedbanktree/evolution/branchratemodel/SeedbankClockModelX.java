package seedbanktree.evolution.branchratemodel;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import seedbanktree.evolution.tree.SeedbankNodeX;

public class SeedbankClockModelX extends BranchRateModel.Base {
	
	final public Input<RealParameter> activeRateParamInput =
            new Input<>("active rate", "the rate parameter associated with active branchs"); 
	
	final public Input<RealParameter> dormantRateParamInput =
            new Input<>("dormant rate", "the rate parameter associated with dormant branchs"); 
	
	final public Input<Tree> treeInput =
            new Input<>("tree", "the tree this relaxed clock is associated with.", Input.Validate.REQUIRED);
	
	private Function meanRate;
	private RealParameter activeRateParameter;
	private RealParameter dormantRateParameter;
	

	@Override
	public void initAndValidate() {
		
		activeRateParameter = activeRateParamInput.get();
		if (activeRateParameter.lowerValueInput.get() == null || activeRateParameter.lowerValueInput.get() < 0.0) {
			activeRateParameter.setLower(0.0);
        }
        if (activeRateParameter.upperValueInput.get() == null || activeRateParameter.upperValueInput.get() < 0.0) {
        	activeRateParameter.setUpper(Double.MAX_VALUE);
        }
        
        dormantRateParameter = dormantRateParamInput.get();
        if (dormantRateParameter.lowerValueInput.get() == null || dormantRateParameter.lowerValueInput.get() < 0.0) {
        	dormantRateParameter.setLower(0.0);
        }
        if (dormantRateParameter.upperValueInput.get() == null || dormantRateParameter.upperValueInput.get() < 0.0) {
        	dormantRateParameter.setUpper(Double.MAX_VALUE);
        }
		
		meanRate = meanRateInput.get();
        if (meanRate == null) {
            meanRate = new RealParameter("1.0");
        }
		
	}
	
	
	@Override
	public double getRateForBranch(Node node) {
		if (((SeedbankNodeX) node).getNodeType() == 0) {
			return dormantRateParameter.getArrayValue() * meanRate.getArrayValue();
		} else if (((SeedbankNodeX) node).getNodeType() == 1) {
			return activeRateParameter.getArrayValue() * meanRate.getArrayValue();
		} else {
			// should not fall through
			throw new RuntimeException("SeedbankNode found with type that isn't 0/1.");
		}
	}
	
}
