package seedbanktree.operators;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import seedbanktree.evolution.tree.SeedbankNode;

public class UniformNode extends SeedbankTreeOperator {

    public Input<Boolean> includeRootInput = new Input<>("includeRoot",
            "Allow modification of root node.", false);
    
    public Input<Double> rootScaleFactorInput = new Input<>("rootScaleFactor",
            "Root scale factor.", 0.9);

    /**
     * Change the node height and return the hastings ratio.
     *
     * @return log of Hastings Ratio
     */
    @Override
    public double proposal() {
        // Randomly select event on tree:
        int event = Randomizer.nextInt(sbTree.getInternalNodeCount() + sbTree.getTotalNumberOfChanges());
        
        
        SeedbankNode node = null;
        int changeIdx = -1;
        if (event<sbTree.getInternalNodeCount()) {
            node = (SeedbankNode)sbTree.getNode(sbTree.getLeafNodeCount() + event);
            if (!includeRootInput.get() && node.isRoot())
                return Double.NEGATIVE_INFINITY;
        } else {
            event -= sbTree.getInternalNodeCount();
            for (Node thisNode : sbTree.getNodesAsArray()) {
                if (thisNode.isRoot())
                    continue;                
                if (event<((SeedbankNode)thisNode).getChangeCount()) {
                    node = (SeedbankNode)thisNode;
                    changeIdx = event;
                    break;
                }
                event -= ((SeedbankNode)thisNode).getChangeCount();
            }
        }
        
        if (node == null)
            throw new IllegalStateException("Event selection loop fell through!");
        
        if (changeIdx==-1) {
            if (node.isRoot()) {
                // Scale distance root and closest event
                double tmin = Math.max(((SeedbankNode)node.getLeft()).getFinalChangeTime(),
                        ((SeedbankNode)node.getRight()).getFinalChangeTime());
                
                double u = Randomizer.nextDouble();
                double f = u*rootScaleFactorInput.get()
                        + (1-u)/rootScaleFactorInput.get();
                
                double tnew = tmin + f*(node.getHeight()-tmin);
                
                node.setHeight(tnew);
                return -Math.log(f);
            } else {
                // Reposition node randomly between closest events
                double tmin = Math.max(((SeedbankNode)node.getLeft()).getFinalChangeTime(),
                        ((SeedbankNode)node.getRight()).getFinalChangeTime());
                
                double tmax = node.getChangeCount()>0
                        ? node.getChangeTime(0)
                        : node.getParent().getHeight();
                
                double u = Randomizer.nextDouble();
                double tnew = u*tmin + (1.0-u)*tmax;
                
                node.setHeight(tnew);
                return 0.0;
            }
        } else {
            double tmin, tmax;
            if (changeIdx+1<node.getChangeCount()) 
                tmax = node.getChangeTime(changeIdx+1);
            else
                tmax = node.getParent().getHeight();
            
            if (changeIdx-1<0)
                tmin = node.getHeight();
            else
                tmin = node.getChangeTime(changeIdx-1);
            
            double u = Randomizer.nextDouble();
            double tnew = u*tmin + (1-u)*tmax;
            
            node.setChangeTime(changeIdx, tnew);
            return 0.0;
        }
        
    }

}
