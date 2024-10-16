package seedbanktree.operators;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import seedbanktree.evolution.tree.SeedbankTree;

public class LambdaScaler extends UniformizationRetypeOperator {

    @Override
    public void initAndValidate() {
    	
    	super.initAndValidate();

        final int indsDim = InputUtil.get(indicatorsInput, this).getDimension();
        final int lambdasDim = InputUtil.get(lambdasInput, this).getDimension();
        if (!(indsDim == lambdasDim)) {
            throw new IllegalArgumentException("indicators dimension not compatible with lambdas dimension");
        }
    }

    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {

        try {

            final RealParameter lambdas = (RealParameter) InputUtil.get(lambdasInput, this);
            final SeedbankTree sbTree = (SeedbankTree) InputUtil.get(seedbankTreeInput, this);
            assert lambdas.getLower() != null && lambdas.getUpper() != null;

            // which position to scale
            final int index;
            final IntegerParameter indicators = (IntegerParameter) InputUtil.get(indicatorsInput, this);
            final int dimCount = indicators.getDimension();
            final Integer[] indicator = indicators.getValues();

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
                return Double.NEGATIVE_INFINITY; // no active indicators
            }
            
            double old_logp = getBranchTypeProb(sbTree.getNode(index));
            double retype_logp = retypeBranch(sbTree.getNode(index));
            if (retype_logp == Double.NEGATIVE_INFINITY)
            	return Double.NEGATIVE_INFINITY;
            
            return old_logp - retype_logp;

        } catch (Exception e) {
            // whatever went wrong, we want to abort this operation...
            return Double.NEGATIVE_INFINITY;
        }
    }
}