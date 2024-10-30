package seedbanktree.operators;

import java.util.Arrays;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import org.jblas.MatrixFunctions;

import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import seedbanktree.evolution.tree.SeedbankNode;
import seedbanktree.evolution.tree.TransitionModel;

public abstract class UniformizationRetypeOperator extends SeedbankTreeOperator {
    
    public UniformizationRetypeOperator() {
    	transitionModelInput.setRule(Validate.REQUIRED);
    }
    
    /**
     * Exception used to signal non-existence of allowed type sequence
     * between node types.
     */
    protected class NoValidPathException extends Exception {
        @Override
        public String getMessage() {
            return "No valid valid type sequence exists between chosen nodes.";
        }
    }
    
    /**
     * Sample the number of virtual events to occur along branch.
     * 
     * General strategy here is to:
     * 1. Draw u from Unif(0,1),
     * 2. Starting from zero, evaluate P(n leq 0|a,b) up until n=thresh
     * or P(n leq 0|a,b)>u.
     * 3. If P(n leq 0|a,b) has exceeded u, use that n. If not, use rejection
     * sampling to draw conditional on n being >= thresh.
     * 
     * @param typeStart Type at start (bottom) of branch
     * @param typeEnd Type at end (top) of branch
     * @param muL Expected unconditioned number of virtual events
     * @param Pba Probability of final type given start type
     * @param transitionModel Transition model to use.
     * @return number of virtual events.
     */
    private int drawEventCount(int typeStart, int typeEnd, double muL, double Pba,
            TransitionModel transitionModel) {
    	
        int nVirt = 0;

        double u = Randomizer.nextDouble();
        double P_low_given_ab = 0.0; // a
        double acc = - muL - Math.log(Pba); 
        double log_muL = Math.log(muL);
        
        do {
        	P_low_given_ab += Math.exp(Math.log(transitionModel.getRpowN(nVirt).get(typeStart, typeEnd)) + acc);
        	
            if (P_low_given_ab>u)
                return nVirt;
            
            nVirt += 1;
            acc += log_muL - Math.log(nVirt);

        } while (transitionModel.RpowSteadyN()<0 || nVirt<transitionModel.RpowSteadyN());
        
        int thresh = nVirt;
        do {
            nVirt = (int) Randomizer.nextPoisson(muL);
        } while (nVirt < thresh);

        return nVirt;
    }
    
    /**
     * Retype branch between srcNode and its parent. Uses the combined
     * uniformization/forward-backward approach of Fearnhead and Sherlock (2006)
     * to condition on both the beginning and end states.
     *
     * @param srcNode
     * @return Probability of new state.
     * @throws seedbanktree.operators.UniformizationRetypeOperator.NoValidPathException
     */
    protected double retypeBranch(Node srcNode) throws NoValidPathException {
        
        Node srcNodeP = srcNode.getParent();
        double t_srcNode = srcNode.getHeight();
        double t_srcNodeP = srcNodeP.getHeight();

        double L = t_srcNodeP-t_srcNode;

        int type_srcNode = ((SeedbankNode)srcNode).getNodeType();
        int type_srcNodeP = ((SeedbankNode)srcNodeP).getNodeType();

        // Pre-calculate some stuff:
        double muL = transitionModel.getMu()*L;
        
        double Pba = MatrixFunctions.expm(
        		transitionModel.getQ()
                .mul(L)).get(type_srcNode,type_srcNodeP);
        
        // Abort if transition is impossible.
        if (Pba == 0.0)
            throw new NoValidPathException();
        
        // Catch for numerical errors
        if (Pba>1.0 || Pba<0.0) {
            System.err.println("Warning: matrix exponentiation resulted in rubbish.  Aborting move.");
            return Double.NEGATIVE_INFINITY;
        }
        
        // Select number of virtual events:
        int nVirt = drawEventCount(type_srcNode, type_srcNodeP, muL, Pba, transitionModel);
        
        if (nVirt<0)
            return Double.NEGATIVE_INFINITY;
        
        // Select times of virtual events:
        double[] times = new double[nVirt];
        for (int i = 0; i<nVirt; i++)
            times[i] = Randomizer.nextDouble()*L+t_srcNode;
        Arrays.sort(times);

        // Sample type changes along branch using FB algorithm:
        int[] types = new int[nVirt];
        int prevType = type_srcNode;
        for (int i = 1; i<=nVirt; i++) {
            
            double u2 = Randomizer.nextDouble()
                    *transitionModel.getRpowN(nVirt-i+1).get(prevType, type_srcNodeP);
            int c;
            boolean fellThrough = true;
            for (c = 0; c<2; c++) {
                u2 -= transitionModel.getR().get(prevType,c)
                        *transitionModel.getRpowN(nVirt-i).get(c,type_srcNodeP);
                if (u2<0.0) {
                    fellThrough = false;
                    break;
                }
            }
            
            // Check for FB algorithm error:
            if (fellThrough) {
                
                double sum1 = transitionModel.getRpowN(nVirt-i+1).get(prevType, type_srcNodeP);
                double sum2 = transitionModel.getR().get(prevType, 0)
                        		* transitionModel.getRpowN(nVirt-i).get(0, type_srcNodeP);
                sum2 += transitionModel.getR().get(prevType, 1)
                        * transitionModel.getRpowN(nVirt-i).get(1, type_srcNodeP);
                
                System.err.println("Warning: FB algorithm failure.  Aborting move."
                		+ " sum1: " + sum1 + ", sum2: " + sum2);
                return Double.NEGATIVE_INFINITY;
            }
            types[i-1] = c;
            prevType = c;
        }

        double logProb = 0.0;
        
        // Add non-virtual type changes to branch, calculating probability
        // of path conditional on start type:
        ((SeedbankNode)srcNode).clearChanges();
        prevType = type_srcNode;
        double prevTime = t_srcNode;
        for (int i = 0; i<nVirt; i++) {

            if (types[i] != prevType) {
                // Add change to branch:
                ((SeedbankNode)srcNode).addChange(types[i], times[i]);

                // Add probability contribution:
                logProb += transitionModel.getQ().get(prevType, prevType)*(times[i]-prevTime)
                        +Math.log(transitionModel.getQ().get(prevType, types[i]));

                prevType = types[i];
                prevTime = times[i];
            }
        }
        logProb += transitionModel.getQ().get(prevType, prevType)*(t_srcNodeP-prevTime);
        
        // Adjust probability to account for end condition:
        logProb -= Math.log(Pba);

        // Return probability of path given boundary conditions:
        return logProb;
    }
    
    /**
     * Obtain probability of the current type-change path above srcNode.
     *
     * @param srcNode
     * @return Path probability.
     */
    protected double getBranchTypeProb(Node srcNode) {
        
    	TransitionModel transitionModel = transitionModelInput.get();

        double logProb = 0.0;

        Node srcNodeP = srcNode.getParent();
        double t_srcNode = srcNode.getHeight();
        double t_srcNodeP = srcNodeP.getHeight();
        double L = t_srcNodeP-t_srcNode;
        int col_srcNode = ((SeedbankNode)srcNode).getNodeType();
        int col_srcNodeP = ((SeedbankNode)srcNodeP).getNodeType();

        // Probability of branch conditional on start type:
        double lastTime = t_srcNode;
        int lastCol = col_srcNode;
        for (int i = 0; i<((SeedbankNode)srcNode).getChangeCount(); i++) {
            double thisTime = ((SeedbankNode)srcNode).getChangeTime(i);
            int thisCol = ((SeedbankNode)srcNode).getChangeType(i);

            logProb += (thisTime-lastTime)*transitionModel.getQ().get(lastCol, lastCol)
                    +Math.log(transitionModel.getQ().get(lastCol, thisCol));

            lastTime = thisTime;
            lastCol = thisCol;
        }
        logProb += (t_srcNodeP-lastTime)*transitionModel.getQ().get(lastCol, lastCol);

        // Adjust to account for end condition of path:
        double Pba = MatrixFunctions.expm(
        		transitionModel.getQ().mul(L)).get(col_srcNode, col_srcNodeP);
        
        // Catch for numerical errors:
        if (Pba>1.0 || Pba < 0.0) {
            System.err.println("Warning: matrix exponentiation resulted in rubbish.  Aborting move.");
            return Double.NEGATIVE_INFINITY;
        }
        
        logProb -= Math.log(Pba);
                
        return logProb;
    }
    
}
