package seedbanktree.evolution.tree;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.BEASTInterface;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;

public class TransitionModel extends CalculationNode {
	
	// Inputs
	// seedbank intensity / rate / c
	public Input<Function> transitionRateInput = new Input<> (
			"transitionRate", 
			"Transition rate (c).", 
			Validate.REQUIRED); 

	// relative seedbank population size K
	// Active pop = K * Dormant pop
	public Input<Function> relativeSeedbankSizeInput = new Input<> (
			"relativeSeedbankSize", 
			"Relative seedbank population size, or the ratio of active to dormant population (K).", 
			Validate.REQUIRED); 
	
	// Active population size
    public Input<Function> activePopSizeInput = new Input<>(
            "activePopSize",
            "Active population size (N).",
            Validate.REQUIRED);
	
    //Shadow inputs
    protected Function transitionRate;
    protected Function relativeSeedbankSize;
    protected Function activePopSize;
	
    @Override
    public void initAndValidate() {
    	
    	transitionRate = transitionRateInput.get();
    	relativeSeedbankSize = relativeSeedbankSizeInput.get();
    	activePopSize = activePopSizeInput.get();

        if (transitionRate instanceof RealParameter) {
        	((RealParameter)transitionRate).setLower(Math.max(((RealParameter)transitionRate).getLower(), 0.0));
            ((RealParameter)transitionRate).setUpper(Math.min(((RealParameter)transitionRate).getUpper(), 1.0));
        }
            
        if (relativeSeedbankSize instanceof RealParameter)
        	((RealParameter)relativeSeedbankSize).setLower(Math.max(((RealParameter)relativeSeedbankSize).getLower(), 0.0));

        if (activePopSize instanceof RealParameter)
            ((RealParameter)activePopSize).setLower(Math.max(((RealParameter)activePopSize).getLower(), 0.0));
        
    }
    
    /**
     * Returns rate of migration
     *
     * @param i, from type
     * @param j, to type
     * @return rate
     */
    public double getBackwardRate(int i, int j) {
        if ( i == j ) {
            return 0;
        } else if ( i == 1 && j == 0) {
        	return transitionRate.getArrayValue();
        } else { //( i == 0 && j == 1 )
        	return transitionRate.getArrayValue() * relativeSeedbankSize.getArrayValue();
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
        	return activePopSize.getArrayValue();
        } else { // i == 0
        	return activePopSize.getArrayValue() / relativeSeedbankSize.getArrayValue();
        }
    }
    
    /**
     * @return a list of the unique identifiers for the parameters describing this transition model
     */
    public List<String> getParameterIds() {
    	// TODO: implement this for SeedbankTreeDensity getConditions();
    	// Code referenced beast.base.evolution.tree.coalescent.ConstantPopulation
    	List<String> ids = new ArrayList<>();
    	
    	if (transitionRateInput.get() instanceof BEASTInterface)
            ids.add(((BEASTInterface)transitionRateInput.get()).getID());
    	
    	if (relativeSeedbankSizeInput.get() instanceof BEASTInterface)
            ids.add(((BEASTInterface)relativeSeedbankSizeInput.get()).getID());
    	
    	if (activePopSizeInput.get() instanceof BEASTInterface)
            ids.add(((BEASTInterface)activePopSizeInput.get()).getID());
    	
        return ids;
    }
}
