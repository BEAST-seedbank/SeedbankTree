package seedbanktree.evolution.branchratemodel;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import seedbanktree.evolution.tree.SeedbankNode;

public class SeedbankClockModel extends BranchRateModel.Base {
	
	final public Input<RealParameter> activeRateParamInput =
            new Input<>("active rate", "the rate parameter associated with active branchs"); 
	
	final public Input<RealParameter> dormantRateParamInput =
            new Input<>("dormant rate", "the rate parameter associated with dormant branchs"); 
	
	final public Input<Tree> treeInput =
            new Input<>("tree", "the tree this relaxed clock is associated with.", Input.Validate.REQUIRED);
	
	private Function meanRate;
	private RealParameter activeRateParameter;
	private RealParameter dormantRateParameter;
	private Tree tree;
	
	private boolean recompute = true;
    private double[] rates; //the output rates

	@Override
	public void initAndValidate() {
		
		tree = treeInput.get();;
		
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
	
	private void calcalateRates(Node node) {
		int nodeNumber = node.getNr();
		
		if (node.isRoot()) {
			rates[nodeNumber] = 1;
		} else {
			SeedbankNode sbNode = (SeedbankNode) node;
			double branchLength = node.getParent().getHeight() - node.getHeight();
			double dormantBranchLength = 0;
			double activeBranchLength = 0;
			
			int type = sbNode.getNodeType();
			
			double previousHeight = node.getHeight();
			double changeHeight = 0;
			for (int i=0; i < sbNode.getChangeCount(); i++) {
				changeHeight = sbNode.getChangeTime(i);
				
				if (type == 0) {
					dormantBranchLength += changeHeight - previousHeight;
				} else if (type == 1) {
					activeBranchLength += changeHeight - previousHeight;
				} else {
					// Should not fall through
					throw new RuntimeException("SeedbankNode found with type that isn't 0/1.");
				}
				
				previousHeight = changeHeight;
				type = 1 - type;
			}
			
			if (type == 0) {
				dormantBranchLength += node.getParent().getHeight() - previousHeight;
			} else if (type == 1) {
				activeBranchLength += node.getParent().getHeight() - previousHeight;
			} else {
				// Should not fall through
				throw new RuntimeException("SeedbankNode found with type that isn't 0/1.");
			}
			
			assert (dormantBranchLength + activeBranchLength == branchLength);
			
			double newRate = 0;
			newRate += (dormantBranchLength / branchLength) * 1;
			newRate += (activeBranchLength / branchLength) * 1;
			rates[nodeNumber] = newRate;
			
		}
		
		if (!node.isLeaf()) {
			calcalateRates(node.getLeft());
			calcalateRates(node.getRight());
		}
	}
	
	@Override
	public double getRateForBranch(Node node) {
		synchronized (this) {
    		if (recompute) {
    			calcalateRates(tree.getRoot());
                recompute = false;
			}
        }

        return rates[node.getNr()] * meanRate.getArrayValue();
	}
	
	@Override
    protected boolean requiresRecalculation() {
        // this is only called if any of its inputs is dirty, hence we need to recompute
        recompute = true;
        return true;
    }

    @Override
    protected void store() {
        recompute = true;
        super.store();
    }

    @Override
    protected void restore() {
        recompute = true;
        super.restore();
    }

}
