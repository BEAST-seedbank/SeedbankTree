package seedbanktree.operators;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;

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
        if (destNode.getHeight()>srcNodeParent.getHeight()
                || srcNode.getHeight()>destNodeParent.getHeight())
            return Double.NEGATIVE_INFINITY;
        
        // Record probability of old colours:
        logHR += getBranchTypeProb(srcNode) + getBranchTypeProb(destNode);
//        manualLog("singledormancy/manualLogging.txt", "" + srcNode.getNr());
//        manualLog("singledormancy/manualLogging.txt", "" + destNode.getNr());
        // Make changes to tree topology:
        replace(srcNodeParent, srcNode, destNode);
        replace(destNodeParent, destNode, srcNode);
        
        logHR -= retypeBranch(srcNode) + retypeBranch(destNode);
        
        return logHR;
    }
    
}
