package seedbanktree.evolution.branchratemodel;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.RealParameter;
import seedbanktree.evolution.tree.SeedbankNode;

public class SeedbankClockModelOld extends BranchRateModel.Base {

	final public Input<Function> dormantScalingInput =
			new Input<>("dormantScaling", "the scaling factor applied to the substitution model of dormant branches.", Validate.REQUIRED);
	
	final public Input<Tree> treeInput =
            new Input<>("tree", "the tree this relaxed clock is associated with.", Validate.REQUIRED);
	
	private Function dormantScaling;
	private Tree tree;
	
	private boolean recompute = true;
    private double[] rates; //the output rates

	@Override
	public void initAndValidate() {
		
		dormantScaling = dormantScalingInput.get();
		if (dormantScaling instanceof RealParameter) {
			RealParameter dS = (RealParameter) dormantScaling;
			
			// scaling parameter range [0, 1]
			dS.setBounds(Math.max(0.0, dS.getLower()), Math.min(1.0, dS.getUpper()));
		}
		
		tree = treeInput.get();
		rates = new double[tree.getNodeCount()];
	}
	
	public void calculateRates(Node node) {
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
			newRate += (dormantBranchLength / branchLength) * dormantScaling.getArrayValue();
			newRate += (activeBranchLength / branchLength);
			rates[nodeNumber] = newRate;
		}
		
		if (!node.isLeaf()) {
			calculateRates(node.getLeft());
			calculateRates(node.getRight());
		}
	}
	
	@Override
	public double getRateForBranch(Node node) {
		synchronized (this) {
    		if (recompute) {
    			calculateRates(tree.getRoot());
                recompute = false;
			}
        }
        return rates[node.getNr()];
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
