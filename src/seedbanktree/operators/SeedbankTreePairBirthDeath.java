package seedbanktree.operators;

import beast.base.evolution.tree.Node;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import seedbanktree.evolution.tree.SeedbankNode;
import seedbanktree.evolution.tree.SeedbankTree;

public class SeedbankTreePairBirthDeath extends SeedbankTreeOperator {
	
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
        
        // Complete either pair birth or pair death proposal:
        if (Randomizer.nextDouble()<0.5)
            return birthProposal(selectedNode, edgeNum, n, m);
        else
            return deathProposal(selectedNode, edgeNum, n, m);
	}
	
    /**
     * Type change pair birth proposal.
     * 
     * @param node Node above which selected edge lies
     * @param edgeNum Number of selected edge
     * @param n Number of nodes on tree.
     * @param m Number of type changes currently on tree.
     * @return log of Hastings factor of move.
     */
    private double birthProposal(Node node, int edgeNum, int n, int m) {
        
        SeedbankNode sbNode = (SeedbankNode)node;
        
        int ridx = edgeNum;
        int sidx = edgeNum-1;
        
        double ts, tr;
        int oldEdgeType;
        if (sidx<0) {
            ts = node.getHeight();
            oldEdgeType = sbNode.getNodeType();
        } else {
            ts = sbNode.getChangeTime(sidx);
            oldEdgeType = sbNode.getChangeType(sidx);
        }

        if (ridx>sbNode.getChangeCount()-1)
            tr = node.getParent().getHeight();
        else
            tr = sbNode.getChangeTime(ridx);

        int newEdgeType = 1 - oldEdgeType;
        
        double tau1 = Randomizer.nextDouble()*(tr-ts) + ts;
        double tau2 = Randomizer.nextDouble()*(tr-ts) + ts;
        double tauMin = Math.min(tau1, tau2);
        double tauMax = Math.max(tau1, tau2);
        
        sbNode.insertChange(edgeNum, oldEdgeType, tauMax);
        sbNode.insertChange(edgeNum, newEdgeType, tauMin);
        
//        return Math.log((migModel.getNTypes()-1)*(m + 2*n - 2)*(tr-ts)*(tr-ts))
//                - Math.log(2*(m + 2*n));
        return Math.log((m + 2*n - 2)*(tr-ts)*(tr-ts))
                - Math.log(2*(m + 2*n));
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
        
        if (sidx<-1 || ridx > sbNode.getChangeCount()) // Non-migration
            return Double.NEGATIVE_INFINITY;
        
        double ts, tr;
        int is, ir;
        if (sidx<0) {
            ts = node.getHeight();
            is = sbNode.getNodeType();
        } else {
            ts = sbNode.getChangeTime(sidx);
            is = sbNode.getChangeType(sidx);
        }
        
        if (ridx>sbNode.getChangeCount()-1)
            tr = node.getParent().getHeight();
        else
            tr = sbNode.getChangeTime(ridx);
        ir = sbNode.getChangeType(ridx-1);
        
//        if (is != ir)
//            return Double.NEGATIVE_INFINITY;
        
        sbNode.removeChange(idx);
        sbNode.removeChange(idx);
        
        return Math.log(2*(m + 2*n - 2))
                - Math.log((m+2*n-4)*(tr-ts)*(tr-ts)); //Not sure where this equation is from
    }
	

}
