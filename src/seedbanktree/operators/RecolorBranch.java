package seedbanktree.operators;

import beast.base.core.Description;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import seedbanktree.evolution.tree.SeedbankNode;
import seedbanktree.evolution.tree.SeedbankTree;

@Description("Proposes a new coloring for a randomly selected branch that "
		+ "already has dormancy.")
public class RecolorBranch extends UniformizationRetypeOperator {

//    @Override
//    public void initAndValidate() {
//    	
//    	super.initAndValidate();
//
//        final int etasDim = InputUtil.get(etasInput, this).getDimension();
//        final int lambdasDim = InputUtil.get(lambdasInput, this).getDimension();
//        if (!(etasDim == lambdasDim)) {
//            throw new IllegalArgumentException("indicators dimension not compatible with lambdas dimension");
//        }
//    }
	
    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
	@Override
    public double proposal() {
		Node node;
		do {
            node = sbTree.getNode(Randomizer.nextInt(sbTree.getNodeCount()));
        } while (node.isRoot());
		
		
		double logHR = 0.0;
		logHR += getBranchTypeProb(node);
		try {
			logHR -= retypeBranch(node);
		} catch (NoValidPathException e) {
			return Double.NEGATIVE_INFINITY;
		}
		
		return logHR;
    }

    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    public double proposal2() {

        try {

            final RealParameter lambdas = (RealParameter) InputUtil.get(lambdasInput, this);
            final SeedbankTree sbTree = (SeedbankTree) InputUtil.get(seedbankTreeInput, this);
            assert lambdas.getLower() != null && lambdas.getUpper() != null;

            // which position to scale
            final int index;
            final IntegerParameter etas = (IntegerParameter) InputUtil.get(etasInput, this);
            final int dimCount = etas.getDimension();
            final Integer[] indicator = etas.getValues();

            // available bit locations. there can be hundreds of them. scan list only once.
            final int[] loc = new int[dimCount + 1];
            int locIndex = 0;

            for (int i = 0; i < dimCount; i++) {
                if (indicator[i] == 1) {
                    loc[locIndex] = i;
                    ++locIndex;
                }
            }

            if (locIndex > 0) {
                final int rand = Randomizer.nextInt(locIndex);
                index = loc[rand];
            } else {
            	// no active etas
                return Double.NEGATIVE_INFINITY;
            }
            
            SeedbankNode sbNode = (SeedbankNode) sbTree.getNode(index);
            double old_logp = getBranchTypeProb(sbNode);
            double retype_logp = retypeBranch(sbNode);
            int tries = 0;
            do {
            	retype_logp = retypeBranch(sbNode);
            	tries++;
                if (retype_logp == Double.NEGATIVE_INFINITY || tries >= 10000000) {
                	manualLog("validation/manualLogging.txt", "TRIES: " + tries + "\n");
                	return Double.NEGATIVE_INFINITY;
                }
            } while (sbNode.getChangeCount() == 0);
            
            if (retype_logp == Double.NEGATIVE_INFINITY)
            	return Double.NEGATIVE_INFINITY;
            
            return old_logp - retype_logp;

        } catch (Exception e) {
            // whatever went wrong, we want to abort this operation...
            return Double.NEGATIVE_INFINITY;
        }
    }
}