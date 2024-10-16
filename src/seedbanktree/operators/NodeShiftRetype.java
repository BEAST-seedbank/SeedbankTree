package seedbanktree.operators;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import seedbanktree.evolution.tree.SeedbankNode;

public class NodeShiftRetype extends UniformizationRetypeOperator {
    
    public Input<Boolean> rootOnlyInput = new Input<>("rootOnly",
            "Always select root node for height adjustment.", false);
    
    public Input<Boolean> noRootInput = new Input<>("noRoot",
            "Never select root node for height adjustment.", false);
    
    public Input<Double> rootScaleFactorInput = new Input<>("rootScaleFactor",
            "Scale factor used in root height proposals. (Default 0.8)", 0.8);
    
    @Override
    public void initAndValidate() {
        super.initAndValidate();
        
        if (rootOnlyInput.get() && noRootInput.get())
            throw new IllegalArgumentException("rootOnly and noRoot inputs "
                    + "cannot both be set to true simultaneously.");
    }
    
    @Override
    public double proposal() {
        // Select internal node to adjust:
        Node node;
        if (rootOnlyInput.get())
            node = sbTree.getRoot();
        else
            do {
                node = sbTree.getNode(sbTree.getLeafNodeCount()
                        + Randomizer.nextInt(sbTree.getInternalNodeCount()));
            } while (noRootInput.get() && node.isRoot());
        
//        manualLog("singledormancy/manualLogging.txt", "" + node.getNr());
        // Generate relevant proposal:
        if (node.isRoot())
            return rootProposal(node);
        else
            return nonRootProposal(node);
    }
    
    /**
     * Root node proposal.
     * @param root
     * @return log of HR
     */
    private double rootProposal(Node root) {
        
        double logHR = 0.0;
        
        // Record probability of current typing:
        logHR += getBranchTypeProb(root.getLeft())
                + getBranchTypeProb(root.getRight());
        
        // Select new root height:
        double u = Randomizer.nextDouble();
        double f = u*rootScaleFactorInput.get() + (1-u)/rootScaleFactorInput.get();
        double oldestChildHeight = Math.max(
                root.getLeft().getHeight(),
                root.getRight().getHeight());
        root.setHeight(oldestChildHeight + f*(root.getHeight()-oldestChildHeight));
        logHR -= Math.log(f);
        
        logHR -= retypeBranch(root.getLeft())
		        + retypeBranch(root.getRight());
        
        recalculateLambda(root.getLeft());
        recalculateLambda(root.getRight());
        
        return logHR;
    }
    
    /**
     * Non-root internal node proposal.
     * @param node
     * @return log of HR of move
     */
    private double nonRootProposal(Node node) {
        
        double logHR = 0.0;
        
        // Record probability of current colouring:
        logHR += getBranchTypeProb(node)
                + getBranchTypeProb(node.getLeft())
                + getBranchTypeProb(node.getRight());
        
        // Select new node height:        
        double upperBound = node.getParent().getHeight();
        double lowerBound = Math.max(
                node.getLeft().getHeight(),
                node.getRight().getHeight());
        node.setHeight(lowerBound+(upperBound-lowerBound)*Randomizer.nextDouble());
        
        logHR -= retypeBranch(node)
		        + retypeBranch(node.getLeft())
		        + retypeBranch(node.getRight());
        
        return logHR;        
    }
    
}
