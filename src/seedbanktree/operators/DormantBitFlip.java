package seedbanktree.operators;

import beast.base.core.Description;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import seedbanktree.evolution.tree.SeedbankNode;

@Description("Flips the value of one eta indicator, and then proposes a "
		+ "new coloring for the corresponding branch if eta is now one, "
		+ "or remove all coloring if eta is now zero.")
public class DormantBitFlip extends UniformizationRetypeOperator {
	
	@Override
	public double proposal() {

		final RealParameter lambdas = ((RealParameter) InputUtil.get(lambdasInput, this));
        final IntegerParameter etas = (IntegerParameter) InputUtil.get(etasInput, this);

        final int dim = etas.getDimension();

        double sum = 0.0;
        for (int i = 0; i < dim; i++) {
        	sum += etas.getValue(i);
        }
        
        final int pos = Randomizer.nextInt(dim);
        final int value = etas.getValue(pos);
        
        double logq = 0.0;
        if (value == 0) {
        	//manualLog("singledormancy/manualLogging.txt", "TURN ON " + pos);
            etas.setValue(pos, 1);
            logq = -Math.log((dim - sum) / (sum + 1));
            
            final double l = Randomizer.nextDouble();
            lambdas.setValue(pos, l);
            
            SeedbankNode sbNode = (SeedbankNode)sbTree.getNode(pos);
//            sbNode = (SeedbankNode) ((SeedbankTree)InputUtil.get(seedbankTreeInput, this)).getNode(pos);
            
//            logq += getBranchTypeProb(sbNode);
//            int tries = 0;
//            double retype_logp;
//            do {
//            	retype_logp = retypeBranch(sbNode);
//            	tries++;
//                if (retype_logp == Double.NEGATIVE_INFINITY || tries >= 10000000) {
////                	manualLog("validation/manualLogging.txt", "TRIES: " + tries + "\n");
//                	return Double.NEGATIVE_INFINITY;
//                }
//            } while (sbNode.getChangeCount() == 0);
////            manualLog("validation/manualLogging.txt", "TRIES: " + tries + "\n");
            
            double retype_logp;
			try {
				retype_logp = retypeBranch(sbNode);
			} catch (NoValidPathException e) {
				return Double.NEGATIVE_INFINITY;
			}
			
            if (retype_logp == Double.NEGATIVE_INFINITY)
            	return Double.NEGATIVE_INFINITY;
            logq -= retype_logp;
            
        } else {
        	//manualLog("singledormancy/manualLogging.txt", "TURN OFF " + pos);
        	SeedbankNode sbNode = (SeedbankNode)sbTree.getNode(pos);
        	
        	if (sbNode.getNodeType() == 0) {
        		return Double.NEGATIVE_INFINITY;
        	}
        	
        	logq += getBranchTypeProb(sbNode);
        	
            etas.setValue(pos, 0);
            lambdas.setValue(pos, 0.0);
            sbNode.clearChanges();
            logq = -Math.log(sum / (dim - sum + 1));
        }
        return logq;
	}

}
