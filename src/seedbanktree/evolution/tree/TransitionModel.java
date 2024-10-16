package seedbanktree.evolution.tree;

import java.util.ArrayList;
import java.util.List;

import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

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
	
    // Shadow inputs
    protected Function transitionRate;
    protected Function relativeSeedbankSize;
    protected Function activePopSize;
    
    protected double mu;
    protected DoubleMatrix Q, R;
    protected List<DoubleMatrix> RpowN;
    protected boolean RpowSteady;
        
    // Flag to indicate whether EV decompositions need updating.
    protected boolean dirty;
    
	
    @Override
    public void initAndValidate() {
    	RpowN = new ArrayList<>();
    	
    	transitionRate = transitionRateInput.get();
    	relativeSeedbankSize = relativeSeedbankSizeInput.get();
    	activePopSize = activePopSizeInput.get();

        if (transitionRate instanceof RealParameter) {
        	((RealParameter)transitionRate).setLower(Math.max(((RealParameter)transitionRate).getLower(), 0.0));
        }
            
        if (relativeSeedbankSize instanceof RealParameter)
        	((RealParameter)relativeSeedbankSize).setLower(Math.max(((RealParameter)relativeSeedbankSize).getLower(), 0.0));

        if (activePopSize instanceof RealParameter)
            ((RealParameter)activePopSize).setLower(Math.max(((RealParameter)activePopSize).getLower(), 0.0));
        
        dirty = true;
        updateMatrices();
    }
    
    /**
     * Ensure all local fields including matrices and eigenvalue decomposition
     * objects are consistent with current values held by inputs.
     */
    public void updateMatrices()  {

        if (!dirty)
            return;

        mu = 0.0;
        Q = new DoubleMatrix(2, 2);
        
        Q.put(0, 1, getBackwardRate(0, 1));
        Q.put(0, 0, -getBackwardRate(0, 1));
        Q.put(1, 0, getBackwardRate(1, 0));
        Q.put(1, 1, -getBackwardRate(1, 0));
        
        mu = Math.max(Q.get(0, 1), Q.get(1, 0));


        // Set up uniformized backward transition rate matrix:
        R = Q.mul(1.0/mu).add(DoubleMatrix.eye(2));
        
        // Clear cached powers of R and steady state flag:
        RpowN.clear();
        RpowSteady = false;
        
        // Power sequences initially contain R^0 = I
        RpowN.add(DoubleMatrix.eye(2));

        dirty = false;
    }
    
    /**
     * Power above which R is known to be steady.
     * 
     * @param symmetric
     * @return index of first known steady element.
     */
    public int RpowSteadyN() {
        if (RpowSteady)
            return RpowN.size();
        else
            return -1;
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
    
    public double getK() {
    	return relativeSeedbankSize.getArrayValue();
    }
    

    public double getMu() {
        updateMatrices();
        return mu;
    }
    

    public DoubleMatrix getR() {
        updateMatrices();
        return R;
    }
    

    public DoubleMatrix getQ() {
        updateMatrices();
        return Q;
    }
    
    public DoubleMatrix getRpowN(int n) {
        updateMatrices();
        
        List <DoubleMatrix> matPowerList;
        DoubleMatrix mat;
        
        matPowerList = RpowN;
        mat = R;
        
        if (n>=matPowerList.size()) {
                
            // Steady state of matrix iteration already reached
            if (RpowSteady) {
                return matPowerList.get(matPowerList.size()-1);
            }
                
            int startN = matPowerList.size();
            for (int i=startN; i<=n; i++) {
                matPowerList.add(matPowerList.get(i-1).mmul(mat));
                                    
                // Occasionally check whether matrix iteration has reached steady state
                if (i%10 == 0) {
                    double maxDiff = 0.0;
                    for (double el : matPowerList.get(i).sub(matPowerList.get(i-1)).toArray())
                        maxDiff = Math.max(maxDiff, Math.abs(el));
                        
                    if (!(maxDiff>0)) {
                        RpowSteady = true;                    
                        return matPowerList.get(i);
                    }
                }
            }
        }
        return matPowerList.get(n);
    }

    /*
     * CalculationNode implementations.
     */
    
    @Override
    protected boolean requiresRecalculation() {
        // we only get here if something is dirty
        dirty = true;
        return true;
    }

    @Override
    protected void restore() {
        dirty = true;
        super.restore();
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
    
    /**
     * Main for debugging.
     *
     * @param args
     */
    public static void main (String [] args) {
        
        int n=10;
        DoubleMatrix Q = new DoubleMatrix(n, n);
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                Q.put(i, j, i*n+j);
            }
        }
        MatrixFunctions.expm(Q.mul(0.001)).print();
        Q.print();
        
    }
}
