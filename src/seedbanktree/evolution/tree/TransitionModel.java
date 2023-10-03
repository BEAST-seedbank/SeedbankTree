package seedbanktree.evolution.tree;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;

public class TransitionModel extends CalculationNode {
	// seedbank intensity c / rate
	public Input<Function> rateInput = new Input<> ("rate", "Transition rate."); 

	// relative active population size K
	public Input<Function> KInput = new Input<> ("K", "Ratio of active to dormant population."); 
	
    public Input<Function> activeSizeInput = new Input<>(
            "activeSize",
            "Active population size.",
            Validate.REQUIRED);
	
	Function rate;
//	Double relativeActiveSize;
	Function K;
	Function activeSize;
	
	public TransitionModel() { }
	
    @Override
    public void initAndValidate() {
    	
    	rate = rateInput.get();
//    	relativeActiveSize = relativeActiveSizeInput.get();
    	K = KInput.get();
    	activeSize = activeSizeInput.get();

        if (rate instanceof RealParameter)
            ((RealParameter)rate).setLower(Math.max(((RealParameter)rate).getLower(), 0.0));
        
        if (K instanceof RealParameter)
        ((RealParameter)K).setLower(Math.max(((RealParameter)K).getLower(), 0.0));

        if (activeSize instanceof RealParameter)
            ((RealParameter)activeSize).setLower(Math.max(((RealParameter)activeSize).getLower(), 0.0));
    }
    
    /**
     * Obtain element of rate matrix for migration model for use in likelihood
     * calculation.  (May be switched to zero in BSSVS calculation.)
     *
     * @param i
     * @param j
     * @return Rate matrix element.
     */
    public double getBackwardRate(int i, int j) {
        if (i==j) {
            return 0;
        } else if (i==1 && j==0) {
        	return rate.getArrayValue();
        } else { //(i==0 && j==1)
        	return rate.getArrayValue() * K.getArrayValue();
        } 
    }
    
    /**
     * Obtain effective population size of particular type/deme.
     *
     * @param i deme index
     * @return Effective population size.
     */
    public double getPopSize(int i) {
        if (i == 1) {
        	return activeSize.getArrayValue();
        } else { // i == 0
        	return activeSize.getArrayValue() / K.getArrayValue();
        }
    }
}
