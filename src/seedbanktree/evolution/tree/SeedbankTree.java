package seedbanktree.evolution.tree;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

import com.google.common.collect.Lists;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.DiscreteStatistics;
import beast.base.util.Randomizer;
import beastfx.app.inputeditor.BeautiDoc;


public class SeedbankTree extends Tree {
	
	// Inputs
	public Input<String> typeLabelInput = new Input<>("typeLabel", "Label for type traits (default 'type')", "type");
	
    public Input<TraitSet> typeTraitInput = new Input<>("typeTrait", "Type trait set.  Used only by BEAUti.");
    
    public Input<String> activeTypeNameInput = new Input<>(
            "activeTypeName", "Name of active type.", "active", Validate.OPTIONAL);
    
    public Input<String> dormantTypeNameInput = new Input<>(
            "activeTypeName", "Name of active type.", "active", Validate.OPTIONAL);

    
	// Shadow inputs
	protected String typeLabel;
	protected TraitSet typeTraitSet;
	protected String activeTypeName;
	protected String dormantTypeName;
	
	@Override
	public void initAndValidate() {
		// If an initial tree is given as input
		if (m_initial.get() != null && !(this instanceof StateNodeInitialiser)) {
            
            if (!(m_initial.get() instanceof SeedbankTree)) {
                throw new IllegalArgumentException("Attempted to initialise "
                        + "seedbank tree with regular tree object.");
            }
            
            SeedbankTree other = (SeedbankTree) m_initial.get();
            root = other.root.copy();
            nodeCount = other.nodeCount;
            internalNodeCount = other.internalNodeCount;
            leafNodeCount = other.leafNodeCount;
        }
		
        if (nodeCount < 0) {
            if (m_taxonset.get() != null) {
                makeCaterpillar(0, 1, false);
            } else {
                // make dummy tree with a single root node
                root = new SeedbankNode();
                root.setNr(0);
                root.setTree(this);
                nodeCount = 1;
                internalNodeCount = 0;
                leafNodeCount = 1;
            }
        }
        
        // nodeCount could be < 0 if m.taxonset.get() != null but is of size 0
		if (nodeCount >= 0) {
            initArrays();
        }
		
        typeLabel = typeLabelInput.get();
        activeTypeName = activeTypeNameInput.get();
        dormantTypeName = dormantTypeNameInput.get();
        
        processTraits(m_traitList.get());
        
     // Ensure tree is compatible with traits.
        if (hasDateTrait())
            adjustTreeNodeHeights(root);
	}
	
	@Override
	public void makeCaterpillar(final double minInternalHeight, final double step, final boolean finalize) {
        // make a caterpillar
        final List<String> taxa = m_taxonset.get().asStringList();
        Node left = new SeedbankNode();
        left.setNr(0);
        left.setHeight(0);
        left.setID(taxa.get(0));
        for (int i = 1; i < taxa.size(); i++) {
            final Node right = new SeedbankNode();
            right.setNr(i);
            right.setHeight(0);
            right.setID(taxa.get(i));
            final Node parent = new SeedbankNode();
            parent.setNr(taxa.size() + i - 1);
            parent.setHeight(minInternalHeight + i * step);
            left.setParent( parent);
            parent.setLeft(left);
            right.setParent(parent);
            parent.setRight(right);
            left = parent;
        }
        root = left;
        leafNodeCount = taxa.size();
        nodeCount = leafNodeCount * 2 - 1;
        internalNodeCount = leafNodeCount - 1;

        if (finalize) {
            initArrays();
        }
    }
	
    @Override
    protected void processTraits(List<TraitSet> traitList) {
        super.processTraits(traitList);
        
        // Only one way to get typeTraitSet at the moment
        // Record trait set associated with leaf types.
        for (TraitSet traitSet : traitList) {
            if (traitSet.getTraitName().equals(typeLabel)) {
                typeTraitSet = traitSet;
                break;
            }
        }
        	
        // typeTraitSet must exist
        if (typeTraitSet == null) {
        	String errorString = String.format(
        			"A type trait set (with name '%s' and types '%s' and '%s') must be provided", 
        			typeLabel, activeTypeName, dormantTypeName);
        	throw new IllegalArgumentException(errorString); 
        }
    }
	
    /**
     * @return TraitSet with same name as typeLabel.
     */
    public TraitSet getTypeTrait() {
        if (!traitsProcessed)
            processTraits(m_traitList.get());
        
        return typeTraitSet;
    }
    
    /**
     * @return true if TraitSet with same name as typeLabel exists.
     */
    public boolean hasTypeTrait() {
        return getTypeTrait() != null;
    }
    
    /**
     * Specifically set the type trait set for this tree. A null value simply
     * removes the existing trait set.
     *
     * @param traitSet
     */
    public void setTypeTrait(TraitSet traitSet) {
        if (hasTypeTrait()) {
            m_traitList.get().remove(typeTraitSet);
        }

        if (traitSet != null) {
            //m_traitList.setValue(traitSet, this);
            typeTraitInput.setValue(traitSet, this);
        }

        typeTraitSet = traitSet;
    }
    
    /**
     * @return type label to be used in logging.
     */
    public String getTypeLabel() {
        return typeLabel;
    }
    
    /**
     * @param typeName name of type
     * @return numerical index representing type
     */
    public int getTypeIndex(String typeName) {
    	if (typeName == activeTypeName) {
    		return 1;
    	} else if (typeName == dormantTypeName) {
    		return 0;
    	} else 
            throw new IllegalArgumentException("TypeSet does not contain type with name " + typeName);
    }
    
    /**
     * @param typeIdx numerical index representing type
     * @return name of type
     */
    public String getTypeName(int typeIdx) {
    	if (typeIdx == 1)
    		return activeTypeName;
    	else if (typeIdx == 0)
			return dormantTypeName;
		else 
            throw new IllegalArgumentException("typeIdx should be either 0 or 1 and not " + typeIdx);
    }
    
    @Override
    public final void initArrays() {
        // initialise tree-as-array representation + its stored variant
        m_nodes = new SeedbankNode[nodeCount];
        listNodes((SeedbankNode)root, (SeedbankNode[])m_nodes);
        m_storedNodes = new SeedbankNode[nodeCount];
        Node copy = root.copy();
        listNodes((SeedbankNode)copy, (SeedbankNode[])m_storedNodes);
    }
    
    /**
     * Convert seedbank tree to array representation.
     *
     * @param node Root of sub-tree to convert.
     * @param nodes Array to populate with tree nodes.
     */
    private void listNodes(SeedbankNode node, SeedbankNode[] nodes) {
        nodes[node.getNr()] = node;
        node.setTree(this);
        if (!node.isLeaf()) {
            listNodes((SeedbankNode)node.getLeft(), nodes);
            if (node.getRight()!=null)
                listNodes((SeedbankNode)node.getRight(), nodes);
        }
    }
    
    /**
     * Deep copy, returns a completely new seedbank tree.
     *
     * @return a deep copy of this seedbank tree
     */
    @Override
    public SeedbankTree copy() {
    	SeedbankTree tree = new SeedbankTree();
        tree.ID = ID;
        tree.index = index;
        tree.root = root.copy();
        tree.nodeCount = nodeCount;
        tree.internalNodeCount = internalNodeCount;
        tree.leafNodeCount = leafNodeCount;
        tree.typeLabel = typeLabel;
        return tree;
    }
    
    /**
     * Copy all values from an existing seedbank tree.
     *
     * @param other
     */
    @Override
    public void assignFrom(StateNode other) {
        SeedbankTree sbTree = (SeedbankTree) other;

        SeedbankNode[] sbNodes = new SeedbankNode[sbTree.getNodeCount()];
        for (int i=0; i<sbTree.getNodeCount(); i++)
        	sbNodes[i] = new SeedbankNode();

        ID = sbTree.ID;
        root = sbNodes[sbTree.root.getNr()];
        root.assignFrom(sbNodes, sbTree.root);
        root.setParent(null);

        nodeCount = sbTree.nodeCount;
        internalNodeCount = sbTree.internalNodeCount;
        leafNodeCount = sbTree.leafNodeCount;
        initArrays();
    }
    
    /**
     * Copy all values aside from IDs from an existing seedbank tree.
     * 
     * @param other
     */
    @Override
    public void assignFromFragile(StateNode other) {
    	SeedbankTree sbTree = (SeedbankTree) other;
        if (m_nodes == null) {
            initArrays();
        }
        root = m_nodes[sbTree.root.getNr()];
        Node[] otherNodes = sbTree.m_nodes;
        int iRoot = root.getNr();
        assignFromFragileHelper(0, iRoot, otherNodes);
        
        root.setHeight(otherNodes[iRoot].getHeight());
        root.setParent(null);
        
        SeedbankNode sbRoot = (SeedbankNode)root;
        sbRoot.nodeType = ((SeedbankNode)(otherNodes[iRoot])).nodeType;
        sbRoot.changeTimes.clear();
        sbRoot.changeTypes.clear();
        sbRoot.nTypeChanges = 0;
        
        if (otherNodes[iRoot].getLeft() != null) {
            root.setLeft(m_nodes[otherNodes[iRoot].getLeft().getNr()]);
        } else {
            root.setLeft(null);
        }
        if (otherNodes[iRoot].getRight() != null) {
            root.setRight(m_nodes[otherNodes[iRoot].getRight().getNr()]);
        } else {
            root.setRight(null);
        }
        assignFromFragileHelper(iRoot + 1, nodeCount, otherNodes);
    }

    /**
     * helper to assignFromFragile
     */
    private void assignFromFragileHelper(int iStart, int iEnd, Node[] otherNodes) {
        for (int i = iStart; i < iEnd; i++) {
        	SeedbankNode sink = (SeedbankNode)m_nodes[i];
        	SeedbankNode src = (SeedbankNode)otherNodes[i];
            
            sink.setHeight(src.getHeight());
            sink.setParent(m_nodes[src.getParent().getNr()]);
            
            sink.nTypeChanges = src.nTypeChanges;
            sink.changeTimes.clear();
            sink.changeTimes.addAll(src.changeTimes);
            sink.changeTypes.clear();
            sink.changeTypes.addAll(src.changeTypes);
            sink.nodeType = src.nodeType;
            
            if (src.getLeft() != null) {
                sink.setLeft(m_nodes[src.getLeft().getNr()]);
                if (src.getRight() != null) {
                    sink.setRight(m_nodes[src.getRight().getNr()]);
                } else {
                    sink.setRight(null);
                }
            }
        }
    }
    
    /**
     * Check whether typing and timing of tree are sensible.
     * 
     * @return true if types and times are "valid"
     */
    public boolean isValid() {
        return timesAreValid(root) && typesAreValid(root);
    }
    
    private boolean timesAreValid(Node node) {
        for (Node child : node.getChildren()) {
            double lastHeight = node.getHeight();
            for (int idx=((SeedbankNode)child).getChangeCount()-1; idx>=0; idx--) { 
                double thisHeight = ((SeedbankNode)child).getChangeTime(idx);
                if (thisHeight>lastHeight)
                    return false;
                lastHeight = thisHeight;
            }
            if (child.getHeight()>lastHeight)
                return false;
            
            if (!timesAreValid(child))
                return false;
        }
        
        return true;
    }
    
    private boolean typesAreValid(Node node) {
    	assert (node instanceof SeedbankNode);
    	
    	if (!node.isLeaf()) { 
    		if (((SeedbankNode)node).getNodeType() != 1)
    			return false;
    		
    		for (Node child : node.getChildren()) {
                if (((SeedbankNode)child).getFinalType() != 1)
                    return false;
                
                int lastType=1;
                for (int idx=((SeedbankNode)child).getChangeCount()-2; idx>=0; idx--) { 
                    int thisType = ((SeedbankNode)child).getChangeType(idx);
                    if (thisType == lastType) // types have to alternate
                        return false;
                    lastType = thisType;
                }
                
                if (!typesAreValid(child))
                    return false;
            }
    	}
        return true;
    }
    
    /**
     * Obtain total number of type changes along nodes on tree.
     * 
     * @return total change count
     */
    public int getTotalNumberOfChanges() {
        int count = 0;        

        for (Node node : m_nodes) {
            if (node.isRoot())
                continue;
            
            count += ((SeedbankNode)node).getChangeCount();
        }
        
        return count;
    }
}