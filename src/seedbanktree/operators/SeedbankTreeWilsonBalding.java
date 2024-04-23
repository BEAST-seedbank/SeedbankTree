package seedbanktree.operators;

import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import seedbanktree.evolution.tree.SeedbankNode;
import seedbanktree.evolution.tree.SeedbankTree;

public class SeedbankTreeWilsonBalding extends SeedbankTreeOperator {

    @Override
    public double proposal() {
    	final SeedbankTree sbTree = (SeedbankTree)InputUtil.get(seedbankTreeInput, this);

    	if (sbTree.getLeafNodeCount() < 3) {
    		throw new IllegalStateException("Tree is too small.");
    	}

        // get random non-root src node
        final int nodeCount = sbTree.getNodeCount();
        SeedbankNode src;
        do {
        	src = (SeedbankNode) sbTree.getNode(Randomizer.nextInt(nodeCount));
        } while (src.isRoot());
        final SeedbankNode srcP = (SeedbankNode) src.getParent();

        // random dest node to insert i above
        SeedbankNode dest;
        SeedbankNode destP;

        // dest is not root and is higher up
        do {
            dest = (SeedbankNode) sbTree.getNode(Randomizer.nextInt(nodeCount));
            destP = (SeedbankNode) dest.getParent();
        } while ((destP == null || destP.getHeight() <= src.getHeight()) || (src.getNr() == dest.getNr()));

        // disallow moves that change the root.
        if (srcP.isRoot()) {
            return Double.NEGATIVE_INFINITY;
        }

        final int srcPnr = srcP.getNr();
        final int destPnr = destP.getNr();
        if ( destPnr == srcPnr || dest.getNr() == srcPnr || destPnr == src.getNr()) {
//        	System.out.println("OPERATOR FAILURE: NODE SELECTION ERROR");
            return Double.NEGATIVE_INFINITY;
        }

        final SeedbankNode sibling = (SeedbankNode) getOtherChild(srcP, src);

        // collect active range for dest branch
        double newMinAge = Math.max(src.getHeight(), dest.getHeight());
        double newRange = 0;
        double lastTime = dest.getHeight();
        for (int i=0; i < dest.getChangeCount(); i++) {
        	double currTime = dest.getChangeTime(i);
        	if (currTime >= newMinAge && dest.getChangeType(i) == 0) {
        		newRange += currTime - Math.max(lastTime, newMinAge);
        	}
        	lastTime = currTime;
        }
        newRange += destP.getHeight() - Math.max(lastTime, newMinAge);
        
        // collect active range for sibling branch
        double oldMinAge = Math.max(src.getHeight(), sibling.getHeight());
        double oldRange = 0;
        lastTime = sibling.getHeight();
        for (int i=0; i < sibling.getChangeCount(); i++) {
        	double currTime = sibling.getChangeTime(i);
        	if (currTime >= oldMinAge && sibling.getChangeType(i) == 0) {
        		oldRange += currTime - Math.max(lastTime, oldMinAge);
        	}
        	lastTime = currTime;
        }
        oldRange += srcP.getHeight() - Math.max(lastTime, oldMinAge);
        
        if (oldRange == 0 || newRange == 0) {
            // This happens when some branch lengths are zero.
            // If oldRange = 0, hastingsRatio == Double.POSITIVE_INFINITY and
            // node i can be catapulted anywhere in the tree, resulting in
            // very bad trees that are always accepted.
            // For symmetry, newRange = 0 should therefore be ruled out as well
            return Double.NEGATIVE_INFINITY;
        }
        
        double hastingsRatio = newRange / Math.abs(oldRange);
        double newAge = Randomizer.nextDouble() * newRange;
                
        // Find actual new age
        double realNewAge;
        
        lastTime = dest.getHeight();
        for (int i=0; i < dest.getChangeCount(); i++) {
        	double currTime = dest.getChangeTime(i);
        	if (currTime >= newMinAge && dest.getChangeType(i) == 0) {
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
        
    	// 1. remove edges <destP, dest>, <srcP, sibling>, <srcPP, srcP>
        // 2. add edges <destP, srcP>, <srcP, dest>, <srcPP, sibling>

        // disconnect p
        final SeedbankNode srcPP = (SeedbankNode) srcP.getParent();
        
        // Implement tree changes:
        disconnectBranch(src);
        connectBranch(src, dest, realNewAge);

        return Math.log(hastingsRatio);
    }


} // class WilsonBalding
