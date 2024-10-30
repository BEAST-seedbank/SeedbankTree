package seedbanktree.operators;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import seedbanktree.evolution.tree.SeedbankNode;

public class DormantBitFlipOld extends UniformizationRetypeOperator {
	
	@Override
	public double proposal() {

		final RealParameter lambdas = ((RealParameter) InputUtil.get(lambdasInput, this));
        final IntegerParameter indicators = (IntegerParameter) InputUtil.get(etasInput, this);

        final int dim = indicators.getDimension();

        double sum = 0.0;
        for (int i = 0; i < dim; i++) {
            if (indicators.getValue(i) == 1) sum += 1;
        }
        
        final int pos = Randomizer.nextInt(dim);
        final int value = indicators.getValue(pos);
        
        double logq = 0.0;
        if (value == 0) {
            indicators.setValue(pos, 1);
            logq = -Math.log((dim - sum) / (sum + 1));
            
            final double l = Randomizer.nextDouble();
            lambdas.setValue(pos, l);
            
            SeedbankNode sbNode = (SeedbankNode)sbTree.getNode(pos);
            
            double retype_logp = constrainedRetypeBranch(sbNode);
            if (retype_logp == Double.NEGATIVE_INFINITY)
            	return Double.NEGATIVE_INFINITY;
            logq -= retype_logp;
            
        } else {
        	SeedbankNode sbNode = (SeedbankNode)sbTree.getNode(pos);
        	
        	if (sbNode.getNodeType() == 0) {
        		return Double.NEGATIVE_INFINITY;
        	}
        	
            indicators.setValue(pos, 0);
            lambdas.setValue(pos, 0.0);
            sbNode.clearChanges();
            logq = -Math.log(sum / (dim - sum + 1));
        }
        return logq;
	}

}
