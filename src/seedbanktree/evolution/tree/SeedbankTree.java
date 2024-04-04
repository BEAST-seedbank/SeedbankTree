package seedbanktree.evolution.tree;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import com.google.common.collect.Lists;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;


public class SeedbankTree extends Tree {
	
	// Inputs
	public Input<String> typeLabelInput = new Input<>("typeLabel", "Label for type traits (default 'type')", "type");
	
//    public Input<TraitSet> typeTraitInput = new Input<>("typeTrait", "Type trait set.  Used only by BEAUti.");
    
    public Input<String> activeTypeNameInput = new Input<>(
            "activeTypeName", "Name of active type.", "active", Validate.OPTIONAL);
    
    public Input<String> dormantTypeNameInput = new Input<>(
            "dormantTypeName", "Name of dormant type.", "dormant", Validate.OPTIONAL);

    
	// Shadow inputs
	protected String typeLabel;
	protected TraitSet typeTraitSet;
	protected String activeTypeName;
	protected String dormantTypeName;
	
	// Constructors
	// Default constructor for beast use
	public SeedbankTree() {}
	
	// For use in initFromFlatTree()
	public SeedbankTree(Node rootNode) {
        
        if (!(rootNode instanceof SeedbankNode))
            throw new IllegalArgumentException("Attempted to instantiate "
                    + "seedbank tree with regular root node.");
        
        setRoot(rootNode);
        initArrays();
    }
	
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
        if (!traitsProcessed) {
        	try {
        		processTraits(m_traitList.get());
        	} catch (IllegalArgumentException e) {}
        }   
        
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

        m_traitList.get().add(traitSet);
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
    	if (typeName.equals(activeTypeName)) {
    		return 1;
    	} else if (typeName.equals(dormantTypeName)) {
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
    	// check root
    	if (root.getLength() != 0.0 || ((SeedbankNode)root).getChangeCount() != 0)
    		return false;
    	
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
    	SeedbankNode sbNode = (SeedbankNode)node;
    	
    	if (!sbNode.isLeaf()) { 
    		if ((sbNode).getNodeType() != 1)
    			return false;
    		
    		for (Node child : node.getChildren()) {
    			assert (child instanceof SeedbankNode);
    			SeedbankNode sbChild = (SeedbankNode)child;
    			
    			switch (sbChild.getChangeCount()) {
    				case 0:
    					if (sbChild.getNodeType() != 1)
    						return false;
    					break;
    					
    				case 1:
    					if (sbChild.getNodeType() != 0 || sbChild.getChangeType(0) != 1)
    						return false;
    					break;
    					
    				default:
    					if (sbChild.getFinalType() != 1)
    	                    return false;
    	                
    	                int lastType=1;
    	                for (int idx=sbChild.getChangeCount()-2; idx>=0; idx--) { 
    	                    int thisType = sbChild.getChangeType(idx);
    	                    if (thisType == lastType) // types have to alternate
    	                        return false;
    	                    lastType = thisType;
    	                }
    	                
    	                if (lastType == sbChild.getNodeType())
    	                	return false;
    			}
                
                if (!typesAreValid(child))
                    return false;
            }
    	}
    	
        return true;
    }
    
    /**
     * Generates a new seedbank tree in which the colours along the branches are
     * indicated by the traits of single-child nodes.
     *
     * This method is useful for interfacing trees coloured externally using the
     * a ColouredTree instance with methods designed to act on trees coloured
     * using single-child nodes and their metadata fields.
     *
     * Caveat: assumes more than one node exists on tree (i.e. leaf != root)
     *
     * @param useTypeStrings whether to use descriptive type strings
     * @return Flattened tree.
     */
    public Tree getFlattenedTree(boolean useTypeStrings) {

        // Create new tree to modify.  Note that copy() doesn't
        // initialise the node array lists, so initArrays() must
        // be called manually.
        Tree flatTree = copy();
        flatTree.initArrays();

        int nextNodeNr = getNodeCount();
        Node colourChangeNode;

        // Iterate over nodes
        for (Node node : getNodesAsArray()) {
 
            SeedbankNode sbNode = (SeedbankNode)node;

            int nodeNum = node.getNr();

            // Root node has no additional changes
            if (node.isRoot()) {
                flatTree.getNode(nodeNum).setMetaData(typeLabel,
                        ((SeedbankNode)(node.getLeft())).getFinalType());
                continue;
            }

            Node startNode = flatTree.getNode(nodeNum);

            startNode.setMetaData(typeLabel,
                    ((SeedbankNode)node).getNodeType());
            if(useTypeStrings)
                startNode.metaDataString = String.format("%s=\"%s\"",
                    typeLabel, getTypeName(sbNode.getNodeType()));
            else
                startNode.metaDataString = String.format("%s=%d",
                    typeLabel, sbNode.getNodeType());

            Node endNode = startNode.getParent();
            
            endNode.setMetaData(typeLabel,
                    ((SeedbankNode)node.getParent()).getNodeType());
            if(useTypeStrings)
                endNode.metaDataString = String.format("%s=\"%s\"",
                    typeLabel, getTypeName(((SeedbankNode)node.getParent()).getNodeType()));
            else
                endNode.metaDataString = String.format("%s=%d",
                    typeLabel, ((SeedbankNode)node.getParent()).getNodeType());

            Node branchNode = startNode;
            for (int i = 0; i<sbNode.getChangeCount(); i++) {

                // Create and label new node:
                colourChangeNode = new SeedbankNode();
                colourChangeNode.setNr(nextNodeNr);
                colourChangeNode.setID(String.valueOf(nextNodeNr));
                nextNodeNr += 1;

                // Connect to child and parent:
                branchNode.setParent(colourChangeNode);
                colourChangeNode.addChild(branchNode);

                // Ensure height and colour trait are set:
                colourChangeNode.setHeight(sbNode.getChangeTime(i));
                colourChangeNode.setMetaData(typeLabel,
                		sbNode.getChangeType(i));
                if (useTypeStrings)
                    colourChangeNode.metaDataString = String.format("%s=\"%s\"",
                        typeLabel, getTypeName(sbNode.getChangeType(i)));
                else
                    colourChangeNode.metaDataString = String.format("%s=%d",
                        typeLabel, sbNode.getChangeType(i));

                // Update branchNode:
                branchNode = colourChangeNode;
            }

            // Ensure final branchNode is connected to the original parent:
            branchNode.setParent(endNode);
            if (endNode.getLeft()==startNode)
                endNode.setLeft(branchNode);
            else
                endNode.setRight(branchNode);
        }

        return flatTree;
    }

    /**
     * Initialise colours and tree topology from Tree object in which colour
     * changes are marked by single-child nodes and colours are stored in
     * meta-data tags. Node numbers of non-singleton nodes in flat tree
     * are preserved.
     *
     * @param flatTree
     * @param takeNrsFromFlatTree 
     * @throws java.lang.Exception 
     */
    public void initFromFlatTree(Tree flatTree, boolean takeNrsFromFlatTree) {

        // Build new coloured tree:

        List<Node> activeFlatTreeNodes = new ArrayList<>();
        List<Node> nextActiveFlatTreeNodes = new ArrayList<>();
        List<SeedbankNode> activeTreeNodes = new ArrayList<>();
        List<SeedbankNode> nextActiveTreeNodes = new ArrayList<>();

        // Populate active node lists with root:
        activeFlatTreeNodes.add(flatTree.getRoot());
        SeedbankNode newRoot = new SeedbankNode();
        activeTreeNodes.add(newRoot);
        
        // Initialise counter used to number leaves when takeNrsFromFlatTree
        // is false:
        int nextNr = 0;

        while (!activeFlatTreeNodes.isEmpty()) {

            nextActiveFlatTreeNodes.clear();
            nextActiveTreeNodes.clear();

            for (int idx = 0; idx<activeFlatTreeNodes.size(); idx++) {
                Node flatTreeNode = activeFlatTreeNodes.get(idx);
                SeedbankNode treeNode = activeTreeNodes.get(idx);

                List<Integer> colours = new ArrayList<>();
                List<Double> times = new ArrayList<>();

                while (flatTreeNode.getChildCount()==1) {
                    Object typeObject = flatTreeNode.getMetaData(typeLabel);
                    int col;
                    if (typeObject instanceof Integer)
                        col = (int) typeObject;
                    else if (typeObject instanceof Double)
                        col = (int) Math.round((double)typeObject);
                    else if (typeObject instanceof String) {
                        try {
                            col = Integer.parseInt((String) typeObject);
                        } catch (NumberFormatException ex) {
                            col = getTypeIndex((String) typeObject);
                        }
                    } else
                        throw new IllegalArgumentException("Unrecognized type metadata.");
                    colours.add(col);

                    times.add(flatTreeNode.getHeight());

                    flatTreeNode = flatTreeNode.getLeft();
                }

                // Order changes to being from youngest to oldest:
                colours = Lists.reverse(colours);
                times = Lists.reverse(times);

                switch (flatTreeNode.getChildCount()) {
                    case 0:
                        // Leaf at base of branch
                        if (takeNrsFromFlatTree) {
                            treeNode.setNr(flatTreeNode.getNr());
                            treeNode.setID(String.valueOf(flatTreeNode.getID()));
                        } else {
                            treeNode.setNr(nextNr);
                            treeNode.setID(String.valueOf(nextNr));
                            nextNr += 1;
                        }
                        break;

                    case 2:
                        // Non-leaf at base of branch
                        nextActiveFlatTreeNodes.add(flatTreeNode.getLeft());
                        nextActiveFlatTreeNodes.add(flatTreeNode.getRight());

                        SeedbankNode daughter = new SeedbankNode();
                        SeedbankNode son = new SeedbankNode();
                        treeNode.addChild(daughter);
                        treeNode.addChild(son);
                        nextActiveTreeNodes.add(daughter);
                        nextActiveTreeNodes.add(son);

                        break;
                }

                // Add type changes to seedbank tree branch:
                for (int i = 0; i<colours.size(); i++)
                    treeNode.addChange(colours.get(i), times.get(i));

                // Set node type at base of seedbank tree branch:
                Object typeObject = flatTreeNode.getMetaData(typeLabel);
                int nodeType;
                if (typeObject instanceof Integer)
                    nodeType = (int)typeObject;
                else if (typeObject instanceof Double)
                    nodeType = (int)Math.round((Double)typeObject);
                else if (typeObject instanceof String) {
                    try {
                        nodeType = Integer.parseInt((String) typeObject);
                    } catch (NumberFormatException ex) {
                        nodeType = getTypeIndex((String) typeObject);
                    }
                } else
                    throw new IllegalArgumentException("Unrecognised type metadata.");

                treeNode.setNodeType(nodeType);

                // Set node height:
                treeNode.setHeight(flatTreeNode.getHeight());
            }

            // Replace old active node lists with new:
            activeFlatTreeNodes.clear();
            activeFlatTreeNodes.addAll(nextActiveFlatTreeNodes);

            activeTreeNodes.clear();
            activeTreeNodes.addAll(nextActiveTreeNodes);

        }
        
        
        // Number internal nodes:
        numberInternalNodes(newRoot, newRoot.getAllLeafNodes().size());
        
        // Assign tree topology:
        assignFromWithoutID(new SeedbankTree(newRoot));
        initArrays();
        
    }
    
    /**
     * Helper method used by initFromFlattenedTree to assign sensible node numbers
     * to each internal node.  This is a post-order traversal, meaning the
     * root is given the largest number.
     * 
     * @param node
     * @param nextNr
     * @return 
     */
    private int numberInternalNodes(Node node, int nextNr) {
        if (node.isLeaf())
            return nextNr;
        
        for (Node child : node.getChildren())
            nextNr = numberInternalNodes(child, nextNr);
        
        node.setNr(nextNr);
        node.setID(String.valueOf(nextNr));
        
        return nextNr+1;
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
    

    /**
     * Stores all data and node data associated with this tree
     */
    @Override
    protected void store() {
        storedRoot = m_storedNodes[root.getNr()];
        int iRoot = root.getNr();

        storeNodes(0, iRoot);
        
        storedRoot.setHeight( m_nodes[iRoot].getHeight());
        storedRoot.setParent(null);

        if (root.getLeft()!=null)
            storedRoot.setLeft(m_storedNodes[root.getLeft().getNr()]);
        else
            storedRoot.setLeft(null);
        if (root.getRight()!=null)
            storedRoot.setRight(m_storedNodes[root.getRight().getNr()]);
        else
            storedRoot.setRight(null);
        
        SeedbankNode sbStoredRoot = (SeedbankNode)storedRoot;
        sbStoredRoot.changeTimes.clear();
        sbStoredRoot.changeTimes.addAll(((SeedbankNode)m_nodes[iRoot]).changeTimes);

        sbStoredRoot.changeTypes.clear();
        sbStoredRoot.changeTypes.addAll(((SeedbankNode)m_nodes[iRoot]).changeTypes);
        
        sbStoredRoot.nTypeChanges = ((SeedbankNode)m_nodes[iRoot]).nTypeChanges;
        sbStoredRoot.nodeType = ((SeedbankNode)m_nodes[iRoot]).nodeType;
        
        storeNodes(iRoot+1, nodeCount);
    }

    /**
     * helper to store *
     */
    private void storeNodes(int iStart, int iEnd) {
        for (int i = iStart; i<iEnd; i++) {
            SeedbankNode sink = (SeedbankNode)m_storedNodes[i];
            SeedbankNode src = (SeedbankNode)m_nodes[i];
            
            sink.setHeight( src.getHeight());
            sink.setParent(m_storedNodes[src.getParent().getNr()]);
            
            if (src.getLeft()!=null) {
                sink.setLeft(m_storedNodes[src.getLeft().getNr()]);
                if (src.getRight()!=null)
                    sink.setRight(m_storedNodes[src.getRight().getNr()]);
                else
                    sink.setRight(null);
            }
            
            sink.changeTimes.clear();
            sink.changeTimes.addAll(src.changeTimes);
            
            sink.changeTypes.clear();
            sink.changeTypes.addAll(src.changeTypes);
            
            sink.nTypeChanges = src.nTypeChanges;
            sink.nodeType = src.nodeType;
        }
    }
    
    /**
     * Return string representation of seedbank tree.
     * 
     * @return Seedbank tree string in Newick format.
     */
    @Override
    public String toString() {

        // Behaves differently if writing a state file
        StackTraceElement[] ste = Thread.currentThread().getStackTrace();
        if (ste[2].getMethodName().equals("toXML")) {            
            // Use toShortNewick to generate Newick string without taxon labels
            return getFlattenedTree(false).getRoot().toShortNewick(true);
        } else{
            return getFlattenedTree(true).getRoot().toSortedNewick(new int[1], true);
        }
    }
    
    /////////////////////////////////////////////////
    // Methods implementing the Loggable interface //
    /////////////////////////////////////////////////
    @Override
    public void init(PrintStream printStream) {

        printStream.println("#NEXUS\n");
        printStream.println("Begin taxa;");
        printStream.println("\tDimensions ntax="+getLeafNodeCount()+";");
        printStream.println("\t\tTaxlabels");
        for (int i = 0; i<getLeafNodeCount(); i++)
            printStream.println("\t\t\t"+getNodesAsArray()[i].getID());
        printStream.println("\t\t\t;");
        printStream.println("End;");

        printStream.println("Begin trees;");
        printStream.println("\tTranslate");
        for (int i = 0; i<getLeafNodeCount(); i++) {
            printStream.print("\t\t\t"+(getNodesAsArray()[i].getNr()+1)
                    +" "+getNodesAsArray()[i].getID());
            if (i<getLeafNodeCount()-1)
                printStream.print(",");
            printStream.print("\n");
        }
        printStream.print(";");
    }

    @Override
    public void log(long i, PrintStream printStream) {
        printStream.print("tree STATE_"+i+" = ");
        printStream.print(toString());
        printStream.print(";");
    }

    @Override
    public void close(PrintStream printStream) {
        printStream.println("End;");
    }
    
    /**
     * Helper function for testing
     */
    public void describeSelf() {
    	System.out.println("DESCRIBING SELF");
    	System.out.println("I am tree " + this.getID());
    	System.out.println(String.format("TypeLabel %s activeTN %s dormantTN %s", typeLabel, activeTypeName, dormantTypeName));
    	System.out.println("traitsProcessed: " + traitsProcessed + " hasTypeTrait(): " + hasTypeTrait());
    	System.out.println(getTypeTrait().traitsInput.get());
    	System.out.println(String.format("nodeCount %d | internalNodeCount %d | leafCount %d", nodeCount, internalNodeCount, leafNodeCount));
    	System.out.println("Iterating through m_nodes -- length: " + m_nodes.length);
    	for (int i=0; i<m_nodes.length; i++) {
    		SeedbankNode n = (SeedbankNode)m_nodes[i];
    		System.out.println("---------");
    		System.out.println("id: " + n.getID());
    		System.out.println("height: " + n.getHeight());
    		System.out.println("branch length: " + n.getLength());
    		System.out.println("num type changes: " + n.getChangeCount());
    	}

    }
}