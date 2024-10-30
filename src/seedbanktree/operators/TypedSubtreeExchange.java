package seedbanktree.operators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;

@Description("Swaps the locations of two nodes in the seedbank tree and "
		+ "recolors the branches between them and their parents.")
public class TypedSubtreeExchange  extends UniformizationRetypeOperator {
    
    public Input<Boolean> isNarrowInput = new Input<Boolean>("isNarrow",
            "Whether or not to use narrow exchange. (Default true.)", true);

    @Override
    public double proposal() {
        double logHR = 0.0;

        // Select source and destination nodes:
        
        Node srcNode, srcNodeParent, destNode, destNodeParent;
        if (isNarrowInput.get()) {
            
            // Narrow exchange selection:
            do {
                srcNode = sbTree.getNode(Randomizer.nextInt(sbTree.getNodeCount()));
            } while (srcNode.isRoot() || srcNode.getParent().isRoot());
            srcNodeParent = srcNode.getParent();            
            destNode = getOtherChild(srcNodeParent.getParent(), srcNodeParent);
            destNodeParent = destNode.getParent();
            
        } else {
            
            // Wide exchange selection:
            do {
                srcNode = sbTree.getNode(Randomizer.nextInt(sbTree.getNodeCount()));
            } while (srcNode.isRoot());
            srcNodeParent = srcNode.getParent();
            do {
                destNode = sbTree.getNode(Randomizer.nextInt(sbTree.getNodeCount()));
            } while(destNode == srcNode
                    || destNode.isRoot()
                    || destNode.getParent() == srcNode.getParent());
            destNodeParent = destNode.getParent();
        }
        
        // Reject if substitution would result in negative branch lengths:
        if (destNode.getHeight() >= srcNodeParent.getHeight()
                || srcNode.getHeight() >= destNodeParent.getHeight())
            return Double.NEGATIVE_INFINITY;
        
        // Record probability of old colours:
        logHR += getBranchTypeProb(srcNode) + getBranchTypeProb(destNode);

        // Make changes to tree topology:
        replace(srcNodeParent, srcNode, destNode);
        replace(destNodeParent, destNode, srcNode);
        
        try {
			logHR -= retypeBranch(srcNode) + retypeBranch(destNode);
		} catch (NoValidPathException e) {
			return Double.NEGATIVE_INFINITY;
		}
        
        return logHR;
    }
    
}
