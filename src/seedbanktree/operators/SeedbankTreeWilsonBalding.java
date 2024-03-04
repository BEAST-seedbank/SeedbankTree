package seedbanktree.operators;

import org.jblas.MatrixFunctions;

import beast.base.evolution.tree.Node;
import beast.base.util.Randomizer;
import multitypetree.evolution.tree.MigrationModel;
import multitypetree.evolution.tree.MultiTypeNode;
import multitypetree.operators.UniformizationRetypeOperator.NoValidPathException;
import seedbanktree.evolution.tree.SeedbankNode;

public class SeedbankTreeWilsonBalding extends SeedbankTreeOperator {

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
        Node srcNodeP = srcNode.getParent();
        Node srcNodeS = getOtherChild(srcNodeP, srcNode);
        double t_srcNode = srcNode.getHeight();
        double t_srcNodeP = srcNodeP.getHeight();
        double t_srcNodeS = srcNodeS.getHeight();

        // Select destination branch node:
        Node destNode;
        do {
            destNode = sbTree.getNode(Randomizer.nextInt(sbTree.getNodeCount()));
        } while (invalidDestNode(srcNode, destNode));
        Node destNodeP = destNode.getParent();
        double t_destNode = destNode.getHeight();

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
            mtTree.setRoot(srcNodeP);

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
            connectBranch(srcNode, destNode, newTime);
            srcNodeS.setParent(null);
            mtTree.setRoot(srcNodeS);

            // Recolour new branch:
            try {
                logHR -= retypeBranch(srcNode);
            } catch (NoValidPathException e) {
                return Double.NEGATIVE_INFINITY;
            }

            // Return HR:
            logHR += Math.log(t_destNodeP-Math.max(t_srcNode, t_destNode))
                    -Math.log(alpha*t_srcNodeS)
                    -(1.0/alpha)*(oldTime/t_srcNodeS-1.0);

            return logHR;
        }
        
        // NON-ROOT MOVE

        double logHR = 0.0;
        SeedbankNode destNode = (SeedbankNode) destNode;

        // Incorporate probability of current colouring.
        logHR += getBranchTypeProb(srcNode);

        // Record srcNode grandmother height:
        double t_srcNodeG = srcNodeP.getParent().getHeight();

        // Choose height of new attachment point:
//        double min_newTime = Math.max(t_destNode, t_srcNode);
        double t_destNodeP = destNodeP.getHeight();
//        double span = t_destNodeP-min_newTime;
//        double newTime = min_newTime+span*Randomizer.nextDouble();
        
        double span;
        double newTime;
        if (destNode.getChangeCount() == 0) {
        	double min_newTime = Math.max(t_destNode, t_srcNode);
        	span = t_destNodeP - min_newTime;
        	newTime = min_newTime+span*Randomizer.nextDouble();
        } else {
	        double prev = t_destNode;
	        for (int i = 0; i < destNode.getChangeCount(); i++) {
	        	int type = destNode.getChangeType(i);
	        	double t = destNode.getChangeTime(i);
	        	
	        	if (t_srcNode > t) {
	        		prev = t;
	        		continue;
	        	}
	        	
	        	if (type == 0) {
	        		if (t_srcNode <= t) {
	        			span += t - t_srcNode;
	        		} else {
	        			span += t - prev;
	        		}
	        	}
	        	prev = t;
	        }
	        
	        span += t_destNodeP - prev;
	        
	        double reverse = span*(1-Randomizer.nextDouble());
	        
	        prev = t_destNodeP; 
	        for (int i = destNode.getChangeCount()-1; i >= 0; i--) {
	        	int type = destNode.getChangeType(i);
	        	double t = destNode.getChangeTime(i);
	        	
	        	if (type == 1) {
	        		if (prev - t >= reverse) {
	        			newTime = prev - reverse;
	        			break;
	        		}
	        		reverse -= prev - t;
	        	}
	        	prev = t;
	        }
	        
        }        

        // Implement tree changes:
        disconnectBranch(srcNode);
        connectBranch(srcNode, destNode, newTime);

//        // Recolour new branch:
//        try {
//            logHR -= retypeBranch(srcNode);
//        } catch (NoValidPathException e) {
//            return Double.NEGATIVE_INFINITY;
//        }

        // HR contribution of topology and node height changes:
        // TODO: May need to be modified
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
     * Obtain probability of the current migratory path above srcNode.
     *
     * @param srcNode
     * @return Path probability.
     */
    protected double getBranchTypeProb(Node srcNode) {
        double logProb = 0.0;

        Node srcNodeP = srcNode.getParent();
        double t_srcNode = srcNode.getHeight();
        double t_srcNodeP = srcNodeP.getHeight();
        double L = t_srcNodeP-t_srcNode;
        int col_srcNode = ((SeedbankNode)srcNode).getNodeType();
        int col_srcNodeP = ((SeedbankNode)srcNodeP).getNodeType();

        // Probability of branch conditional on start type:
        double lastTime = t_srcNode;
        int lastCol = col_srcNode;
        for (int i = 0; i<((SeedbankNode)srcNode).getChangeCount(); i++) {
            double thisTime = ((SeedbankNode)srcNode).getChangeTime(i);
            int thisCol = ((SeedbankNode)srcNode).getChangeType(i);

            logProb += (thisTime-lastTime)*migrationModel.getQ(sym).get(lastCol, lastCol)
                    +Math.log(migrationModel.getQ(sym).get(lastCol, thisCol));

            lastTime = thisTime;
            lastCol = thisCol;
        }
        logProb += (t_srcNodeP-lastTime)*migrationModel.getQ(sym).get(lastCol, lastCol);

        // Adjust to account for end condition of path:
        double Pba = MatrixFunctions.expm(
                migrationModel.getQ(sym).mul(L)).get(col_srcNode, col_srcNodeP);
        
        // Catch for numerical errors:
        if (Pba>1.0 || Pba < 0.0) {
            System.err.println("Warning: matrix exponentiation resulted in rubbish.  Aborting move.");
            return Double.NEGATIVE_INFINITY;
        }
        
        logProb -= Math.log(Pba);
                
        return logProb;
    }

}
