package seedbanktree.operators;

import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import seedbanktree.evolution.tree.SeedbankNode;

public class TypedWilsonBalding extends UniformizationRetypeOperator {

    @Override
    public double proposal() {
    	
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
        double t_srcNodeS = srcNodeS.getHeight();

        // Select destination branch node:
        SeedbankNode destNode;
        do {
            destNode = (SeedbankNode)sbTree.getNode(Randomizer.nextInt(sbTree.getNodeCount()));
        } while (invalidDestNode(srcNode, destNode));
        Node destNodeP = destNode.getParent();
        double t_destNode = destNode.getHeight();

        
        double logHR = 0.0;

        // Incorporate probability of current colouring.
        logHR += getBranchTypeProb(srcNode);

        // Record srcNode grandmother height:
        double t_srcNodeG = srcNodeP.getParent().getHeight();

        // Choose height of new attachment point:
        double min_newTime = Math.max(t_destNode, t_srcNode);
        double t_destNodeP = destNodeP.getHeight();
        double span = t_destNodeP - min_newTime;
        double newTime = min_newTime + span * Randomizer.nextDouble();

        // Implement tree changes:
        disconnectBranch(srcNode);
        connectBranch(srcNode, destNode, newTime);

        // Recolour new branch AND incorporate probability of new coloring:
        try {
            logHR -= retypeBranch(srcNode);
        } catch (NoValidPathException e) {
            return Double.NEGATIVE_INFINITY;
        }
        
        logHR += Math.log(t_destNodeP-Math.max(t_srcNode, t_destNode))
                -Math.log(t_srcNodeG-Math.max(t_srcNode, t_srcNodeS));
        return logHR;
    }

    /**
     * Returns true if srcNode CANNOT be used for the CWBR move.
     *
     * @param srcNode
     * @return True if srcNode invalid.
     */
    private boolean invalidSrcNode(Node srcNode) {

        if (srcNode.isRoot() || srcNode.getParent().isRoot())
            return true;

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
                || destNode==srcNode.getParent()
                || destNode.getParent()==srcNode.getParent()
                || destNode.isRoot())
            return true;

        Node destNodeP = destNode.getParent();

        if (destNodeP != null && (destNodeP.getHeight() <= srcNode.getHeight()))
            return true;

        return false;
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