package seedbanktree.operators;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import seedbanktree.evolution.tree.SeedbankNode;

public class TypedWilsonBalding extends UniformizationRetypeOperator {

    public Input<Double> alphaInput = new Input<>("alpha",
            "Root height proposal parameter", Validate.REQUIRED);
    private double alpha;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        alpha = alphaInput.get();
    }

    @Override
    public double proposal() {
//    	retypeBranch(sbTree.getNode(0));
//    	if (0 == 0)
//    		return 0.0;
    	
        // Check that operator can be applied to tree:
        if (sbTree.getLeafNodeCount()<3)
            throw new IllegalStateException("Tree too small for"
                    +" TypedWilsonBalding operator.");

        // Select source node:
        Node srcNode;
        do {
            srcNode = sbTree.getNode(Randomizer.nextInt(sbTree.getNodeCount()));
        } while (invalidSrcNode(srcNode));
        SeedbankNode srcNodeP = (SeedbankNode) srcNode.getParent();
        SeedbankNode srcNodeS = (SeedbankNode) getOtherChild(srcNodeP, srcNode);
        double t_srcNode = srcNode.getHeight();
        double t_srcNodeP = srcNodeP.getHeight();
        double t_srcNodeS = srcNodeS.getHeight();

        // Select destination branch node:
        SeedbankNode destNode;
        do {
            destNode = (SeedbankNode)sbTree.getNode(Randomizer.nextInt(sbTree.getNodeCount()));
        } while (invalidDestNode(srcNode, destNode));
        Node destNodeP = destNode.getParent();
        double t_destNode = destNode.getHeight();

        if (destNode.isRoot() || srcNodeP.isRoot()) {
        	return Double.NEGATIVE_INFINITY;
        }
        
        // Handle special cases involving root:
        if (destNode.isRoot()) {
        	// FORWARD ROOT MOVE

            double logHR = 0.0;

            // Record probability of current colouring:
            logHR += getBranchTypeProb(srcNode);
            
            // Record srcNode grandmother height:
            double t_srcNodeG = srcNodeP.getParent().getHeight();

            // Choose new root height:
            double newTime = t_destNode+Randomizer.nextExponential(1.0/(alpha*t_destNode));

            // Implement tree changes:
            disconnectBranch(srcNode);
            connectBranchToRoot(srcNode, destNode, newTime);
            sbTree.setRoot(srcNodeP);

            // Recolour root branches:
            try {
                logHR -= retypeRootBranches(srcNode);
            } catch (NoValidPathException e) {
                return Double.NEGATIVE_INFINITY;
            }

            // Return HR:
            logHR += Math.log(alpha*t_destNode)
                    +(1.0/alpha)*(newTime/t_destNode-1.0)
                    -Math.log(t_srcNodeG-Math.max(t_srcNode, t_srcNodeS));
            return logHR;
        }
        if (srcNodeP.isRoot()) {
            // BACKWARD ROOT MOVE
            
            double logHR = 0.0;

            // Incorporate probability of current colouring:
            logHR += getRootBranchTypeProb(srcNode);

            // Record old srcNode parent height
            double oldTime = t_srcNodeP;

            // Choose height of new attachement point:
            double min_newTime = Math.max(t_srcNode, t_destNode);
            double t_destNodeP = destNodeP.getHeight();
            double span = t_destNodeP-min_newTime;
            double newTime = min_newTime+span*Randomizer.nextDouble();

            // Implement tree changes:
            disconnectBranchFromRoot(srcNode);
            connectBranch(srcNode, destNode,newTime);
            srcNodeS.setParent(null);
            sbTree.setRoot(srcNodeS);

            logHR -= retypeBranch(srcNode);

            // Return HR:
            logHR += Math.log(t_destNodeP-Math.max(t_srcNode, t_destNode))
                    -Math.log(alpha*t_srcNodeS)
                    -(1.0/alpha)*(oldTime/t_srcNodeS-1.0);
            return logHR;
        }
        
        // NON-ROOT MOVE

        double logHR = 0.0;

        // Incorporate probability of current colouring.
        logHR += getBranchTypeProb(srcNode);

        // Record srcNode grandmother height:
//        double t_srcNodeG = srcNodeP.getParent().getHeight();

        // Choose height of new attachment point:
//        double min_newTime = Math.max(t_destNode, t_srcNode);
//        double t_destNodeP = destNodeP.getHeight();
//        double span = t_destNodeP-min_newTime;
//        double newTime = min_newTime+span*Randomizer.nextDouble();

//        manualLog("./singledormancy/manualLogging.txt", "src: " + srcNode.getNr());
//        manualLog("./singledormancy/manualLogging.txt", "dest: " + destNode.getNr());
        
        // collect active range for dest branch
        double newMinAge = Math.max(t_srcNode, t_destNode);
        double newRange = 0;
        double lastTime = destNode.getHeight();
        for (int i=0; i < destNode.getChangeCount(); i++) {
        	double currTime = destNode.getChangeTime(i);
        	if (currTime >= newMinAge && destNode.getChangeType(i) == 0) {
        		newRange += currTime - Math.max(lastTime, newMinAge);
        	}
        	lastTime = currTime;
        }
        newRange += destNodeP.getHeight() - Math.max(lastTime, newMinAge);
//        manualLog("./singledormancy/manualLogging.txt", "newRange: " + newRange);
        
        // collect active range for return, which is sibling branch + parent branch
        double oldMinAge = Math.max(t_srcNode, t_srcNodeS);
        double oldRange = 0;
        lastTime = t_srcNodeS;
        for (int i=0; i < srcNodeS.getChangeCount(); i++) {
        	double currTime = srcNodeS.getChangeTime(i);
        	if (currTime >= oldMinAge && srcNodeS.getChangeType(i) == 0) {
        		oldRange += currTime - Math.max(lastTime, oldMinAge);
        	}
        	lastTime = currTime;
        }
        oldRange += t_srcNodeP - Math.max(lastTime, oldMinAge);
//        manualLog("./singledormancy/manualLogging.txt", "oldRange: " + oldRange);
        
        lastTime = t_srcNodeP;
        for (int i=0; i < srcNodeP.getChangeCount(); i++) {
        	double currTime = srcNodeP.getChangeTime(i);
        	if (srcNodeP.getChangeType(i) == 0) {
        		oldRange += currTime - lastTime;
        	}
        	lastTime = currTime;
        }
        oldRange += srcNodeP.getParent().getHeight() - lastTime;
//        manualLog("./singledormancy/manualLogging.txt", "oldRange: " + oldRange);
        
        
        if (oldRange == 0 || newRange == 0) {
            // This happens when some branch lengths are zero.
            // If oldRange = 0, hastingsRatio == Double.POSITIVE_INFINITY and
            // node i can be catapulted anywhere in the tree, resulting in
            // very bad trees that are always accepted.
            // For symmetry, newRange = 0 should therefore be ruled out as well
            return Double.NEGATIVE_INFINITY;
        }
        
        double newAge = Randomizer.nextDouble() * newRange;
        
     // Find actual new age
        double realNewAge;
//        manualLog("./singledormancy/manualLogging.txt", "newMinAge: " + newMinAge);
//        manualLog("./singledormancy/manualLogging.txt", "newAge: " + newAge);
        
        lastTime = t_destNode;
        for (int i=0; i < destNode.getChangeCount(); i++) {
        	double currTime = destNode.getChangeTime(i);
        	if (currTime >= newMinAge && destNode.getChangeType(i) == 0) {
        		double segLen = currTime - Math.max(lastTime, newMinAge);
        		if (newAge == segLen || newAge == 0.0) {
        			return Double.NEGATIVE_INFINITY;
        		} else if (newAge > segLen) {
        			newAge -= segLen;
        		} else {
        			break;
        		}
        	}
        	lastTime = currTime;
        }
        realNewAge = Math.max(lastTime, newMinAge) + newAge;
        
        
        
//        manualLog("./singledormancy/manualLogging.txt", "realNewAge: " + realNewAge);
        
        // Implement tree changes:
        disconnectBranch(srcNode);
        connectBranch(srcNode, destNode, realNewAge);

//        manualLog("./singledormancy/manualLogging.txt", "logHR: " + logHR);
        logHR -= retypeBranch(srcNode);
//        manualLog("./singledormancy/manualLogging.txt", "logHR: " + logHR);
        // HR contribution of topology and node height changes:
        logHR += Math.log(newRange) - Math.log(oldRange);
//        logHR += Math.log(t_destNodeP-Math.max(t_srcNode, t_destNode))
//                -Math.log(t_srcNodeG-Math.max(t_srcNode, t_srcNodeS));
//        manualLog("./singledormancy/manualLogging.txt", "logHR: " + logHR);
        return logHR;
    }

    /**
     * Returns true if srcNode CANNOT be used for the CWBR move.
     *
     * @param srcNode
     * @return True if srcNode invalid.
     */
    private boolean invalidSrcNode(Node srcNode) {

        if (srcNode.isRoot())
            return true;

        Node parent = srcNode.getParent();

        // This check is important for avoiding situations where it is
        // impossible to choose a valid destNode:
        if (parent.isRoot()) {

            Node sister = getOtherChild(parent, srcNode);

            if (sister.isLeaf())
                return true;

            if (srcNode.getHeight()>=sister.getHeight())
                return true;
        }

        return false;
    }

    /**
     * Returns true if destNode CANNOT be used for the CWBR move in conjunction
     * with srcNode.
     *
     * @param srcNode
     * @param destNode
     * @return True if destNode invalid.
     */
    private boolean invalidDestNode(Node srcNode, Node destNode) {

        if (destNode==srcNode
                ||destNode==srcNode.getParent()
                ||destNode.getParent()==srcNode.getParent())
            return true;

        Node destNodeP = destNode.getParent();

        if (destNodeP!=null&&(destNodeP.getHeight()<=srcNode.getHeight()))
            return true;

        return false;
    }
    
    /**
     * Retype branches with nChanges between srcNode and the root (srcNode's
     * parent) and nChangesSister between the root and srcNode's sister.
     *
     * @param srcNode
     * @return Probability of new state.
     */
    private double retypeRootBranches(Node srcNode) throws NoValidPathException {
        
        double logProb = 0.0;

        Node srcNodeP = srcNode.getParent();
        Node srcNodeS = getOtherChild(srcNodeP, srcNode);

        // Select new root colour:
//        ((SeedbankNode)srcNodeP).setNodeType(Randomizer.nextInt(2));
        
        // Incorporate probability of new root colour:
//        logProb += Math.log(1.0/2);

        // Recolour branches conditional on root type:
        logProb += retypeBranch(srcNode);
        logProb += retypeBranch(srcNodeS);


        // Return probability of new colouring given boundary conditions:
        return logProb;
    }
    
    
    /**
     * Obtain joint probability of typing along branches between srcNode and
     * the root, the sister of srcNode and the root, and the node type of the
     * root.
     *
     * @param srcNode
     * @return
     */
    protected double getRootBranchTypeProb(Node srcNode) {
        
        double logProb = 0.0;

        Node srcNodeP = srcNode.getParent();
        Node srcNodeS = getOtherChild(srcNodeP, srcNode);
        
        // Probability of node type:
        logProb += Math.log(1.0/2);

        // Probability of branch types conditional on node types:
        logProb += getBranchTypeProb(srcNode);
        logProb += getBranchTypeProb(srcNodeS);

        return logProb;
    }
    
}