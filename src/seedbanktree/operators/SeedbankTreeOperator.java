package seedbanktree.operators;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import seedbanktree.evolution.tree.SeedbankTree;
import seedbanktree.evolution.tree.TransitionModel;
import seedbanktree.evolution.tree.SeedbankNode;

@Description("This operator generates proposals for a seedbank tree")
public abstract class SeedbankTreeOperator extends Operator {
	
	    final public Input<SeedbankTree> seedbankTreeInput = 
	    		new Input<>("sbTree", "Seedbank tree on which to operate.", Validate.REQUIRED);
            
        final public Input<TransitionModel> transitionModelInput = 
        		new Input<>("transitionModel", "Transition model for active and dormant types", Validate.OPTIONAL);
        
        protected SeedbankTree sbTree;
        protected TransitionModel transitionModel;
        
        @Override
        public void initAndValidate() {
        	sbTree = seedbankTreeInput.get();
        	transitionModel = transitionModelInput.get();
        }

        /* ***********************************************************************
         * The following two methods are copied verbatim from TreeOperator.
         */
        
        /**
         * Obtain the sister of node "child" having parent "parent".
         * 
         * @param parent the parent
         * @param child  the child that you want the sister of
         * @return the other child of the given parent.
         */
        protected Node getOtherChild(Node parent, Node child) {
            if (parent.getLeft().getNr() == child.getNr()) {
                return parent.getRight();
            } else {
                return parent.getLeft();
            }
        }

        /**
         * Replace node "child" with another node.
         *
         * @param node
         * @param child
         * @param replacement
         */
        public void replace(Node node, Node child, Node replacement) {
        	node.removeChild(child);
        	node.addChild(replacement);
            node.makeDirty(Tree.IS_FILTHY);
            replacement.makeDirty(Tree.IS_FILTHY);
        }
        
        /* **********************************************************************/

        /**
         * Disconnect edge <node,node.getParent()> from the tree by joining 
         * node's sister directly to node's grandmother and adding all colour 
         * changes previously on <node.getParent(),node.getParent().getParent()> 
         * to the new branch.
         *
         * @param node
         */
        public void disconnectBranch(Node node) {
        	
            // Check argument validity:
            SeedbankNode parent = (SeedbankNode) node.getParent();
            if (node.isRoot() || parent.isRoot())
                throw new IllegalArgumentException("Illegal argument to "
                        + "disconnectBranch().");

            SeedbankNode sister = (SeedbankNode) getOtherChild(parent, node);

            // Add colour changes originally attached to parent to those attached
            // to node's sister:
            for (int idx = 0; idx < (parent).getChangeCount(); idx++) {
                int colour = parent.getChangeType(idx);
                double time = parent.getChangeTime(idx);
                sister.addChange(colour, time);
            }

            // Implement topology change.
            replace(parent.getParent(), parent, sister);

            // Clear colour changes from parent:
            parent.clearChanges();
            
            // Ensure BEAST knows to update affected likelihoods:
            parent.makeDirty(Tree.IS_FILTHY);
            sister.makeDirty(Tree.IS_FILTHY);
            node.makeDirty(Tree.IS_FILTHY);
        }

        /**
         * Disconnect node from root, discarding all colouring on <node,root> and
         * <node's sister,root>. The node's sister is made the new root.
         *
         * @param node
         */
        public void disconnectBranchFromRoot(Node node) {

            // Check argument validity:
            if (node.isRoot() || !node.getParent().isRoot())
                throw new IllegalArgumentException("Illegal argument to"
                        + " disconnectBranchFromRoot().");

            // Implement topology change:
            Node parent = node.getParent();
            Node sister = getOtherChild(parent, node);
            sister.setParent(null);
            parent.removeChild(sister);

            // Clear colour changes on new root:
            ((SeedbankNode)sister).clearChanges();

            // Ensure BEAST knows to update affected likelihoods:
            parent.makeDirty(Tree.IS_FILTHY);
            sister.makeDirty(Tree.IS_FILTHY);
            node.makeDirty(Tree.IS_FILTHY);
        }

        /**
         * Creates a new branch between node and a new node at time destTime between
         * destBranchBase and its parent. Colour changes are divided between the two
         * new branches created by the split.
         *
         * @param node
         * @param destBranchBase
         * @param destTime
         */
        public void connectBranch(Node node,
        		Node destBranchBase, double destTime) {

            // Check argument validity:
            if (node.isRoot() || destBranchBase.isRoot())
                throw new IllegalArgumentException("Illegal argument to "
                        + "connectBranch().");

            // Obtain existing parent of node and set new time:
            SeedbankNode parent = (SeedbankNode) node.getParent();
            parent.setHeight(destTime);
            
            // Determine where the split comes in the list of colour changes
            // attached to destBranchBase:
            SeedbankNode sbDestBranchBase = (SeedbankNode)destBranchBase;
            int split;
            for (split = 0; split < sbDestBranchBase.getChangeCount(); split++)
                if (sbDestBranchBase.getChangeTime(split) > destTime)
                    break;

            // Divide colour changes between new branches:
            parent.clearChanges();
            for (int idx = split; idx < sbDestBranchBase.getChangeCount(); idx++)
            	parent.addChange(sbDestBranchBase.getChangeType(idx),
            			sbDestBranchBase.getChangeTime(idx));

            sbDestBranchBase.truncateChanges(split);

            // Set colour at split:
            parent.setNodeType(sbDestBranchBase.getFinalType());

            // Implement topology changes:
            replace(destBranchBase.getParent(), destBranchBase, parent);
            destBranchBase.setParent(parent);

            if (parent.getLeft() == node)
                parent.setRight(destBranchBase);
            else if (parent.getRight() == node)
                parent.setLeft(destBranchBase);

            // Ensure BEAST knows to update affected likelihoods:
            node.makeDirty(Tree.IS_FILTHY);
            parent.makeDirty(Tree.IS_FILTHY);
            destBranchBase.makeDirty(Tree.IS_FILTHY);
        }

        /**
         * Set up node's parent as the new root with a height of destTime, with
         * oldRoot as node's new sister.
         *
         * @param node
         * @param oldRoot
         * @param destTime
         */
        public void connectBranchToRoot(Node node, Node oldRoot, double destTime) {

            // Check argument validity:
            if (node.isRoot() || !oldRoot.isRoot())
                throw new IllegalArgumentException("Illegal argument "
                        + "to connectBranchToRoot().");

            // Obtain existing parent of node and set new time:
            Node newRoot = node.getParent();
            newRoot.setHeight(destTime);

            // Implement topology changes:

            newRoot.setParent(null);

            if (newRoot.getLeft() == node)
                newRoot.setRight(oldRoot);
            else if (newRoot.getRight() == node)
                newRoot.setLeft(oldRoot);

            oldRoot.setParent(newRoot);

            // Ensure BEAST knows to recalculate affected likelihood:
            newRoot.makeDirty(Tree.IS_FILTHY);
            oldRoot.makeDirty(Tree.IS_FILTHY);
            node.makeDirty(Tree.IS_FILTHY);
        }
}
