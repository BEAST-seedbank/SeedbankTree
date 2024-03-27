package seedbanktree.operators;

import beast.base.evolution.tree.Node;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import seedbanktree.evolution.tree.SeedbankNode;
import seedbanktree.evolution.tree.SeedbankTree;

public class SeedbankTreeDormantScale extends SeedbankTreeOperator {

	@Override
	public double proposal() {
		final SeedbankTree sbTree = (SeedbankTree)InputUtil.get(seedbankTreeInput, this);
		
		// Get total number of dormant branches
		int d = 0;
		
		for (Node leaf : sbTree.getExternalNodes()) {
			SeedbankNode sbLeaf = (SeedbankNode) leaf;
			if (sbLeaf.getNodeType() == 1)
				d += sbLeaf.getChangeCount() / 2;
			else
				d += (sbLeaf.getChangeCount() + 1) / 2;
        }
		
		for (Node node : sbTree.getInternalNodes()) {
			SeedbankNode sbNode = (SeedbankNode) node;
			d += sbNode.getChangeCount() / 2;
        }
		
		if (d == 0) return Double.NEGATIVE_INFINITY;
		
		// Get dormant edge
		int edgeNum = Randomizer.nextInt(d);
		
		SeedbankNode selectedNode = null;
        for (Node node : sbTree.getNodesAsArray()) {
        	SeedbankNode sbNode = (SeedbankNode) node;
            if (sbNode.isRoot())
                continue;
            
            int isDormant = sbNode.getNodeType() == 0 ? 1 : 0;
            int dormantCount = (sbNode.getChangeCount() + isDormant) / 2;
            if (edgeNum < dormantCount) {
                selectedNode = sbNode;
                break;
            }
            edgeNum -= dormantCount;
        }
        
        // If the branch extends from a dormant sample, then only one end of the interval is changed
        if (selectedNode.getNodeType() == 0 && edgeNum == 0) {
        	double lower = selectedNode.getHeight();
        	double upper;
			if (selectedNode.getChangeCount() == 1) 
        		upper = selectedNode.getParent().getHeight();
        	else 
        		upper = selectedNode.getChangeTime(1);
			
        	double newTime = lower + Randomizer.nextDouble() * (upper - lower);
        	
        	if (newTime == lower) return Double.NEGATIVE_INFINITY;
        	
        	selectedNode.setChangeTime(0, newTime);
        	
        } else {
        	double lower, upper;
        	
        	int isDormant = selectedNode.getNodeType() == 0 ? 1 : 0;
        	int baseIdx = edgeNum * 2 - isDormant;
        	
        	lower = baseIdx == 0 ? selectedNode.getHeight() : selectedNode.getChangeTime(baseIdx - 1);
        	if (selectedNode.getChangeCount() == baseIdx + 2) 
        		upper = selectedNode.getParent().getHeight();
        	else 
        		upper = selectedNode.getChangeTime(baseIdx + 2);
        	
        	double newBaseTime = lower + Randomizer.nextDouble() * (upper - lower);
        	double newEndTime = lower + Randomizer.nextDouble() * (upper - lower);
        	
        	if (newBaseTime == lower || newEndTime == lower || newBaseTime == newEndTime) 
        		return Double.NEGATIVE_INFINITY;
        	
        	if (newBaseTime > newEndTime) {
        		double temp = newBaseTime;
        		newBaseTime = newEndTime;
        		newEndTime = temp;
        	}
        	
        	selectedNode.setChangeTime(baseIdx, newBaseTime);
        	selectedNode.setChangeTime(baseIdx + 1, newEndTime);
        	
        }
		
		return 0.0;
	}
}
