package seedbanktree.operators;

import beast.base.evolution.tree.Node;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import seedbanktree.evolution.tree.SeedbankNode;
import seedbanktree.evolution.tree.SeedbankTree;

public class SeedbankTreePairDeath extends SeedbankTreeOperator {
	
	@Override
	public double proposal() {
		final SeedbankTree sbTree = (SeedbankTree)InputUtil.get(seedbankTreeInput, this);
		
		int n = sbTree.getLeafNodeCount();
        int m = sbTree.getTotalNumberOfChanges();
        
        // Select sub-edge at random:
        int edgeNum = Randomizer.nextInt(2*n - 2 + m);
        
        // Find edge that sub-edge lies on:
        Node selectedNode = null;
        for (Node node : sbTree.getNodesAsArray()) {
            if (node.isRoot())
                continue;

            if (edgeNum<((SeedbankNode)node).getChangeCount()+1) {
                selectedNode = node;
                break;
            }
            edgeNum -= ((SeedbankNode)node).getChangeCount()+1;
        }
        
        return deathProposal(selectedNode, edgeNum, n, m);
	}
	
    /**
     * Colour change pair death proposal.
     * 
     * @param node Node above which selected edge lies
     * @param edgeNum Number of selected edge
     * @param n Number of nodes on tree
     * @param m Number of colour changes currently on tree
     * @return log of Hastings factor of move.
     */
    private double deathProposal(Node node, int edgeNum, int n, int m) {
        
        SeedbankNode sbNode = (SeedbankNode)node;
        
        int idx = edgeNum-1;
        int sidx = edgeNum-2;
        int ridx = edgeNum+1;
        
        if (sidx<-1 || ridx > sbNode.getChangeCount())
            return Double.NEGATIVE_INFINITY;
        
        double ts, tr;
        if (sidx<0) {
            ts = node.getHeight();
        } else {
            ts = sbNode.getChangeTime(sidx);
        }
        
        if (ridx>sbNode.getChangeCount()-1)
            tr = node.getParent().getHeight();
        else
            tr = sbNode.getChangeTime(ridx);
               
        sbNode.removeChange(idx);
        sbNode.removeChange(idx);
        
        return Math.log(2*(m + 2*n - 2))
                - Math.log((m+2*n-4)*(tr-ts)*(tr-ts));
    }
}
