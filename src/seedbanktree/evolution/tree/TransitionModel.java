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
	public Input<Function> rateInput = new Input<> ("rate", "Transition rate."); 

	// relative active population size K
	// Active pop = K * Dormant pop
	public Input<Function> KInput = new Input<> ("K", "Ratio of active to dormant population."); 
	
	// Active population size
    public Input<Function> activeSizeInput = new Input<>(
            "activeSize",
            "Active population size.",
            Validate.REQUIRED);
	
    //Shadow inputs
    protected Function rate;
    protected Function K;
    protected Function activeSize;
	private String activeTypeName;
    private String dormantTypeName;
	
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
        	return rate.getArrayValue();
        } else { //( i == 0 && j == 1 )
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
    
    /**
     * @return a list of the unique identifiers for the parameters describing this transition model
     */
    public List<String> getParameterIds() {
    	// TODO: implement this for SeedbankTreeDensity getConditions();
    	// Code referenced beast.base.evolution.tree.coalescent.ConstantPopulation
    	List<String> ids = new ArrayList<>();
    	
    	if (rateInput.get() instanceof BEASTInterface)
            ids.add(((BEASTInterface)rateInput.get()).getID());
    	
    	if (KInput.get() instanceof BEASTInterface)
            ids.add(((BEASTInterface)KInput.get()).getID());
    	
    	if (activeSizeInput.get() instanceof BEASTInterface)
            ids.add(((BEASTInterface)activeSizeInput.get()).getID());
    	
        return ids;
    }
}
