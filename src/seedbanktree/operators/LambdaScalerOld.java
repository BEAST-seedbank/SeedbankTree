package seedbanktree.operators;

import java.text.DecimalFormat;

import beast.base.core.Input;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import seedbanktree.evolution.tree.SeedbankTree;

public class LambdaScalerOld extends UniformizationRetypeOperator {

    public final Input<Double> scaleFactorInput = new Input<>("scaleFactor", "scaling factor: range from 0 to 1. Close to zero is very large jumps, close to 1.0 is very small jumps.", 0.75);
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);
    final public Input<Double> scaleUpperLimit = new Input<>("upper", "Upper Limit of scale factor", 1.0 - 1e-8);
    final public Input<Double> scaleLowerLimit = new Input<>("lower", "Lower limit of scale factor", 1e-8);

    /**
     * shadows input *
     */
    private double scaleFactor;

    private double upper, lower;

    @Override
    public void initAndValidate() {
    	
    	super.initAndValidate();

        scaleFactor = scaleFactorInput.get();
        upper = scaleUpperLimit.get();
        lower = scaleLowerLimit.get();

        final int indsDim = InputUtil.get(indicatorsInput, this).getDimension();
        final int lambdasDim = InputUtil.get(lambdasInput, this).getDimension();
        if (!(indsDim == lambdasDim)) {
            throw new IllegalArgumentException("indicators dimension not compatible with lambdas dimension");
        }
    }

    protected boolean outsideBounds(final double value, final RealParameter param) {
        final Double l = param.getLower();
        final Double h = param.getUpper();

        return (value < l || value > h);
    }

    protected double getScaler() {
        return (scaleFactor + (Randomizer.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));
    }

    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {

        try {

            final double scale = getScaler();
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
            //manualLog("singledormancy/manualLogging.txt", "scaling " + index);
            final double oldValue = lambdas.getValue(index);

            if (oldValue == 0) {
                // Error: parameter has value 0 and cannot be scaled
            	//manualLog("singledormancy/manualLogging.txt", "oldValue == 0");
                return Double.NEGATIVE_INFINITY;
            }

            final double newValue = scale * oldValue;

            if (outsideBounds(newValue, lambdas)) {
                // reject out of bounds scales
            	//manualLog("singledormancy/manualLogging.txt", "outsideBounds");
                return Double.NEGATIVE_INFINITY;
            }

            lambdas.setValue(index, newValue);
            
            double retype_logp = constrainedRetypeBranch(sbTree.getNode(index));
            if (retype_logp == Double.NEGATIVE_INFINITY)
            	return Double.NEGATIVE_INFINITY;
            
            return -retype_logp;
            
//            SeedbankNode sbNode = (SeedbankNode)sbTree.getNode(index);
//            final double diff = (newValue - oldValue) * sbNode.getLength();
//            double changeTime0 = sbNode.getChangeTime(0);
//            double changeTime1 = sbNode.getChangeTime(1);
//            double timeBelow = changeTime0 - sbNode.getHeight();
//            double timeAbove = sbNode.getParent().getHeight() - changeTime1;
//            
//            double changeTime0New = changeTime0 - (diff * (timeBelow / (timeBelow + timeAbove)));
//            double changeTime1New = changeTime1 + (diff * (timeAbove / (timeBelow + timeAbove)));
//            
//            sbNode.clearChanges();
//            sbNode.addChange(0, changeTime0New);
//            sbNode.addChange(1, changeTime1New);
            
//            return -Math.log(scale);

        } catch (Exception e) {
            // whatever went wrong, we want to abort this operation...
            return Double.NEGATIVE_INFINITY;
        }
    }


    /**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
        if (optimiseInput.get()) {
            double delta = calcDelta(logAlpha);
            delta += Math.log(1.0 / scaleFactor - 1.0);
            setCoercableParameterValue(1.0 / (Math.exp(delta) + 1.0));
        }
    }

    @Override
    public double getCoercableParameterValue() {
        return scaleFactor;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        scaleFactor = Math.max(Math.min(value, upper), lower);
    }

    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        final double sf = Math.pow(scaleFactor, ratio);

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > 0.40) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }

}
