package seedbanktree.operators;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import seedbanktree.evolution.tree.SeedbankNode;
import seedbanktree.evolution.tree.SeedbankTree;

public class DormancyMover extends Operator{

    public final Input<RealParameter> lambdasInput = new Input<>("lambdas", "the parameter to be scaled", Input.Validate.REQUIRED);
    final public Input<IntegerParameter> indicatorsInput = new Input<>("indicators", "indicates which of the dimension " +
            "of the parameters can be scaled.", Input.Validate.REQUIRED);
    
    final public Input<SeedbankTree> sbTreeInput = new Input<>("sbTree", "Seedbank Tree", Validate.REQUIRED);


    @Override
    public void initAndValidate() {
        final IntegerParameter indicators = indicatorsInput.get();
        final int dataDim = lambdasInput.get().getDimension();
        final int indsDim = indicators.getDimension();
        if (!(indsDim == dataDim)) {
            throw new IllegalArgumentException("indicator dimension not compatible from parameter dimension");
        }
    }

    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {

        final RealParameter lambdas = (RealParameter) InputUtil.get(lambdasInput, this);
        final SeedbankTree sbTree = (SeedbankTree) InputUtil.get(sbTreeInput, this);

        // which position to scale
        final int index;
        final IntegerParameter indicators = indicatorsInput.get();
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

        final double l = lambdas.getValue(index);

        if (l == 0) {
            // Error: parameter has value 0 and cannot be scaled
            return Double.NEGATIVE_INFINITY;
        }

        SeedbankNode sbNode = (SeedbankNode)sbTree.getNode(index);
        sbNode.clearChanges();
        
        final double transitionLoc = Randomizer.nextDouble() * (sbNode.getLength() * (1 - l));
        sbNode.addChange(0, sbNode.getHeight() + transitionLoc);
        sbNode.addChange(1, sbNode.getHeight() + transitionLoc + (sbNode.getLength() * l));

        return 0;
    }

}
