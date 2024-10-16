package seedbanktree.operators;

import java.util.Arrays;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

import org.jblas.MatrixFunctions;

import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import seedbanktree.evolution.tree.SeedbankNode;
import seedbanktree.evolution.tree.TransitionModel;

public abstract class UniformizationRetypeOperator extends SeedbankTreeOperator {
	
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
    
    public UniformizationRetypeOperator() {
    	transitionModelInput.setRule(Validate.REQUIRED);
//    	lambdasInput.setRule(Validate.REQUIRED);
//    	indicatorsInput.setRule(Validate.REQUIRED);
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
     * @param migrationModel Migration model to use.
     * @return number of virtual events.
     */
    private int drawEventCount(int typeStart, int typeEnd, double muL, double Pba,
            TransitionModel transitionModel) {
        
    	// lambda * delta in supplement 
    	// ~ m * delta in Vaughan?
    	// ~ mu * lambda in Rodrigue
    	
    	
        int nVirt = 0;

        double u = Randomizer.nextDouble();
        double P_low_given_ab = 0.0; // a
        
        double acc = - muL - Math.log(Pba); 
        // log[e^(-lambda*delta) * (lambda * delta)^n / n!) / Pba]
        // The division by Pba comes from g = p(b|a, lambda)*U being the threshold,
        // dividing by Pba here makes it so only a direct comparison to U is needed later.
        // -muL represents e^(-lambda*delta)
        
        double log_muL = Math.log(muL);
        
        do {
        	// Step 5
        	P_low_given_ab += Math.exp(Math.log(transitionModel.getRpowN(nVirt).get(typeStart, typeEnd)) + acc);
            // Step 6
            if (P_low_given_ab>u)
                return nVirt;
            
            nVirt += 1;
            acc += log_muL - Math.log(nVirt); // Step 5 prep (log of muL^n, 1/n!)
            
            // Step 4
            // while not at steady state
        } while (transitionModel.RpowSteadyN()<0 || nVirt<transitionModel.RpowSteadyN());
        
        int thresh = nVirt;
        
        // P_n_given_ab constant for n>= thresh: only need
        // to sample P(n|n>=thresh)
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
     * @throws multitypetree.operators.UniformizationRetypeOperator.NoValidPathException
     */
    protected double retypeBranch(Node srcNode) {
//    	Tree t = srcNode.getTree();
//    	srcNode = t.getNode(287);
//    	Node srcNodeP_ = srcNode.getParent();
//        double t_srcNode_ = srcNode.getHeight();
//        double t_srcNodeP_ = srcNodeP_.getHeight();
//        double L_ = t_srcNodeP_-t_srcNode_;
//        manualLog("spikeandslab/manualLogging.txt", "" + L_);
//
//        int type_srcNode_ = 1;
//        int type_srcNodeP_ = 1;
//        manualLog("spikeandslab/manualLogging.txt", "" + type_srcNode_ + " " + type_srcNodeP_);
//
//        // Pre-calculate some stuff:
//        double muL_ = transitionModel.getMu()*L_;
//        manualLog("spikeandslab/manualLogging.txt", "" + muL_);
//        
//        double Pba_ = MatrixFunctions.expm(
//        		transitionModel.getQ()
//                .mul(L_)).get(type_srcNode_,type_srcNodeP_);
//        manualLog("spikeandslab/manualLogging.txt", "" + Pba_);
//    	
//    	double[] nVirt_arr = new double[100];
//    	double[] changes_arr = new double[100];
//    	for (int j = 0; j < 10000; j++) {
//        	int nVirt_i = drawEventCount(type_srcNode_, type_srcNodeP_, muL_, Pba_,
//        			transitionModel);
//            nVirt_arr[nVirt_i]++;
//            
//            int[] types_arr = new int[nVirt_i];
//            int prevType_ = type_srcNode_;
//            for (int i = 1; i<=nVirt_i; i++) {
//                
//                double u2 = Randomizer.nextDouble()
//                        *transitionModel.getRpowN(nVirt_i-i+1).get(prevType_, type_srcNodeP_);
//                int c;
//                boolean fellThrough = true;
//                for (c = 0; c<2; c++) {
//                    u2 -= transitionModel.getR().get(prevType_,c)
//                            *transitionModel.getRpowN(nVirt_i-i).get(c,type_srcNodeP_);
//                    if (u2<0.0) {
//                        fellThrough = false;
//                        break;
//                    }
//                }
//                
//                // Check for FB algorithm error:
//                if (fellThrough) {
//                    
//                    double sum1 = transitionModel.getRpowN(nVirt_i-i+1).get(prevType_, type_srcNodeP_);
//                    double sum2 = 0;
//                    for (c = 0; c<2; c++) {
//                        sum2 += transitionModel.getR().get(prevType_,c)
//                                *transitionModel.getRpowN(nVirt_i-i).get(c,type_srcNodeP_);
//                    }
//                    
//                    System.err.println("Warning: FB algorithm failure.  Aborting move.");
//                    return Double.NEGATIVE_INFINITY;
//                }
//
//                types_arr[i-1] = c;
//                prevType_ = c;
//            }
//            int changes=0;
//            for (int i = 1; i < nVirt_i; i++) {
//            	if (types_arr[i] != types_arr[i-1])
//            		changes++;
//            }
//            if (nVirt_i > 0 && types_arr[0] != t_srcNode_) {
//            	changes++;
//            }
//            changes_arr[changes]++;
//        }
//    	
//    	manualLog("spikeandslab/manualLogging.txt", Arrays.toString(nVirt_arr));
//    	manualLog("spikeandslab/manualLogging.txt", Arrays.toString(changes_arr));
//    	
//    	if (0 == 0)
//    		return Double.NEGATIVE_INFINITY;
        
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
        
        // Catch for numerical errors
        if (Pba>1.0 || Pba<0.0) {
            System.err.println("Warning: matrix exponentiation resulted in rubbish.  Aborting move.");
            return Double.NEGATIVE_INFINITY;
        }
        
        // Select number of virtual events:
        int nVirt = drawEventCount(type_srcNode, type_srcNodeP, muL, Pba, transitionModel);
//        manualLog("validation/manualLogging.txt", "muL, nVirt: " + muL +", " + nVirt);
        
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
        recalculateLambda(srcNode);
        
        // Adjust probability to account for end condition:
        logProb -= Math.log(Pba);

        // Return probability of path given boundary conditions:
        return logProb;
    }
    
    /**
     * Obtain probability of the current migratory path above srcNode.
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
    
    
    /**
     * Retype branch between srcNode and its parent. Uses the combined
     * uniformization/forward-backward approach of Fearnhead and Sherlock (2006)
     * to condition on both the beginning and end states, with an additional
     * constraint imposed by the existing lambda value.
     *
     * @param srcNode
     * @return Probability of new state.
     * @throws multitypetree.operators.UniformizationRetypeOperator.NoValidPathException
     */
    protected double constrainedRetypeBranch(Node srcNode) {
    	//manualLog("singledormancy/manualLogging.txt", "constrainedRetypeBranch " + srcNode.getNr());
    	assert (indicatorsInput.get().getValue(srcNode.getNr()) == 1);
    	
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
        
        // Catch for numerical errors
        if (Pba>1.0 || Pba<0.0) {
        	//manualLog("singledormancy/manualLogging.txt", "aborting move");
            System.err.println("Warning: matrix exponentiation resulted in rubbish.  Aborting move.");
            return Double.NEGATIVE_INFINITY;
        }
        
        // Select number of virtual events:
        int nVirt = drawEventCount(type_srcNode, type_srcNodeP, muL, Pba, transitionModel);
        
        if (nVirt<0) {
        	//manualLog("singledormancy/manualLogging.txt", "nVirt < 0, " + nVirt);
            return Double.NEGATIVE_INFINITY;
        }
        
        //manualLog("singledormancy/manualLogging.txt", "muL = " + muL);	
        //manualLog("singledormancy/manualLogging.txt", "Pba = " + Pba);	
        //manualLog("singledormancy/manualLogging.txt", "nVirt = " + nVirt);	

        // Sample type changes along branch using FB algorithm:
        int[] types = new int[nVirt];
        int prevType = type_srcNode;
        for (int i = 1; i<=nVirt; i++) {
        	if (nVirt == 16) {
            	//manualLog("singledormancy/manualLogging.txt", "i, prevType: " + i + ", " + prevType);
            	//manualLog("singledormancy/manualLogging.txt", "prev to parent " + transitionModel.getRpowN(nVirt-i+1).get(prevType, type_srcNodeP));
            	//manualLog("singledormancy/manualLogging.txt", "prevType to 0 " + transitionModel.getR().get(prevType,0));
            	//manualLog("singledormancy/manualLogging.txt", "0 to parent " + transitionModel.getRpowN(nVirt-i).get(0,type_srcNodeP));
            	//manualLog("singledormancy/manualLogging.txt", "prevType to 0 to parent " + transitionModel.getR().get(prevType,0)*transitionModel.getRpowN(nVirt-i).get(0,type_srcNodeP));
            	//manualLog("singledormancy/manualLogging.txt", "prevType to 1 " + transitionModel.getR().get(prevType,1));
            	//manualLog("singledormancy/manualLogging.txt", "1 to parent " + transitionModel.getRpowN(nVirt-i).get(1,type_srcNodeP));
            	//manualLog("singledormancy/manualLogging.txt", "prevType to 1 to parent " + transitionModel.getR().get(prevType,1)*transitionModel.getRpowN(nVirt-i).get(1,type_srcNodeP));
            	//manualLog("singledormancy/manualLogging.txt", "");
            }
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
                //manualLog("singledormancy/manualLogging.txt", "aborting move 2");
                return Double.NEGATIVE_INFINITY;
            }
            types[i-1] = c;
            prevType = c;
        }

        // Count number of branch segments
        int lastType = type_srcNode;
        int numDormant = 0;
        for (int i = 0; i < nVirt; i++) {
        	if (types[i] == 1 && lastType == 0) {
        		numDormant += 1;
        	}
        	lastType = types[i];
        }
        int numActive = numDormant + type_srcNode;
        
        if (numDormant == 0) {
        	//manualLog("singledormancy/manualLogging.txt", "numDormant == 0");
        	return Double.NEGATIVE_INFINITY;
        }
        //manualLog("singledormancy/manualLogging.txt", "L: " + L + ", lambda: " + InputUtil.get(lambdasInput, this).getArrayValue(srcNode.getNr()));
        //manualLog("singledormancy/manualLogging.txt", "types " + Arrays.toString(types));
        //manualLog("singledormancy/manualLogging.txt", "numDormant/Active: " + numDormant + ", " + numActive);
        // Get dormant segment lengths
        double[] dormantBreaks = new double[numDormant-1];
        for (int i = 0; i<numDormant-1; i++)
        	dormantBreaks[i] = Randomizer.nextDouble();
        Arrays.sort(dormantBreaks);
        //manualLog("singledormancy/manualLogging.txt", "dormantBreaks " + Arrays.toString(dormantBreaks));
        
        double dormant_L = L * InputUtil.get(lambdasInput, this).getArrayValue(srcNode.getNr());
        double lastBreak = 0; 
        double[] dormantLengths = new double[numDormant];
        for (int i = 0; i<numDormant-1; i++) {
        	dormantLengths[i] = dormant_L * dormantBreaks[i] - lastBreak;
        	lastBreak = dormant_L * dormantBreaks[i];
        }
    	dormantLengths[numDormant-1] = dormant_L;
//    	if (numDormant >= 2) 
		dormantLengths[numDormant-1] = dormant_L - lastBreak;
    	
    	//manualLog("singledormancy/manualLogging.txt", "dormantLengths " +Arrays.toString(dormantLengths));
        
    	// Get active segment lengths
    	double[] activeBreaks = new double[numActive-1];
    	for (int i = 0; i<numActive-1; i++)
        	activeBreaks[i] = Randomizer.nextDouble();
    	Arrays.sort(activeBreaks);
    	//manualLog("singledormancy/manualLogging.txt", "activeBreaks " +Arrays.toString(activeBreaks));
    	
    	double active_L = L - dormant_L;
    	lastBreak = 0;
        double[] activeLengths = new double[numActive];
        for (int i = 0; i<numActive-1; i++) {
        	activeLengths[i] = active_L * activeBreaks[i] - lastBreak;
        	lastBreak = active_L * activeBreaks[i];
        }
        activeLengths[numActive-1] = active_L;
//        if (numActive >= 2)
    	activeLengths[numActive-1] = active_L - lastBreak;
    	//manualLog("singledormancy/manualLogging.txt", "activeLengths " +Arrays.toString(activeLengths));
        
    	// Combine active and dormant segment lengths
        double[] lengths = new double[numDormant + numActive];
        for (int i = 0; i<numDormant; i++) {
        	lengths[i*2] = type_srcNode == 0 ? dormantLengths[i] : activeLengths[i];
        	lengths[i*2 + 1] = type_srcNode == 0 ? activeLengths[i] : dormantLengths[i];
        }
        if (type_srcNode == 1)
        	lengths[numDormant+numActive-1] = activeLengths[numActive-1];
        
        assert(DoubleStream.of(lengths).sum() == L);
        //manualLog("singledormancy/manualLogging.txt", "lengths " +Arrays.toString(lengths));
        double logProb = 0.0;
        
        // Add non-virtual type changes to branch, calculating probability
        // of path conditional on start type:
        ((SeedbankNode)srcNode).clearChanges();
        prevType = type_srcNode;
        double prevTime = t_srcNode;
        for (int i = 0; i<numDormant+numActive-1; i++) {
            // Add change to branch:
        	double time = prevTime + lengths[i]; 
        	//manualLog("singledormancy/manualLogging.txt", "AddChange " + (1-prevType) + ", " + time);
            ((SeedbankNode)srcNode).addChange(1 - prevType, time);

            // Add probability contribution:
            logProb += transitionModel.getQ().get(prevType, prevType)*(time-prevTime)
                    +Math.log(transitionModel.getQ().get(prevType, 1 - prevType));

            prevType = 1 - prevType;
            prevTime = time;
        }
        logProb += transitionModel.getQ().get(prevType, prevType)*(t_srcNodeP-prevTime);
        recalculateLambda(srcNode);
        
        // Adjust probability to account for end condition:
        logProb -= Math.log(Pba);

        // Return probability of path given boundary conditions:
        return logProb;
    }

}
