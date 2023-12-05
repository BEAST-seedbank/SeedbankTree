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

import seedbanktree.evolution.tree.SeedbankNodeX.NodeEvent;


public class SeedbankTreeX extends Tree implements StateNodeInitialiser {
	
	// Fields
	
	// Inputs MTT
	public Input<String> typeLabelInput = new Input<>("typeLabel", "Label for type traits (default 'type')", "type");
    public Input<TraitSet> typeTraitInput = new Input<>("typeTrait", "Type trait set.  Used only by BEAUti.");
	
	// Non-inputs MTT
	protected String typeLabel;
    protected TraitSet typeTraitSet;
//    typeSet not necessary -- only two types, pre-defined
//    protected TypeSet typeSet;
    
    // Inputs SCMTT
    public Input<TransitionModel> transitionModelInput = new Input<>(
            "transitionMode,",
            "Transition model to use in simulator.",
            Validate.REQUIRED);
    
    public Input<IntegerParameter> leafTypesInput = new Input<>(
            "leafTypes",
            "Types of leaf nodes.");
    
    public Input<String> outputFileNameInput = new Input<>(
            "outputFileName", "Optional name of file to write simulated "
                    + "tree to.");
    
    // Non-inputs SCMTT
    protected TransitionModel transitionModel;
    
    private List<Integer> leafTypes;
    private List<String> leafNames;
    private List<Double> leafTimes;
    private int nLeaves;
    
    /*
     * Other private fields and classes:
     */
    private abstract class SBEvent {

        double time;
        int fromType, toType;

    }

    private class CoalescenceEvent extends SBEvent {

        public CoalescenceEvent(int type, double time) {
        	// TODO: coalescent events don't need from and to... rework events?
            this.fromType = type;
            this.time = time;
        }
    }

    private class MigrationEvent extends SBEvent {

        public MigrationEvent(int fromType, int toType, double time) {
            this.fromType = fromType;
            this.toType = toType;
            this.time = time;
        }
    }
    
    private class NullEvent extends SBEvent {
        public NullEvent() {
            this.time = Double.POSITIVE_INFINITY;
        }
    }
    
	
	
    public SeedbankTreeX() {}
	
	@Override
	public void initAndValidate() {
		if (m_initial.get() != null) {
            
            if (!(m_initial.get() instanceof SeedbankTreeX)) {
                throw new IllegalArgumentException("Attempted to initialise "
                        + "seedbank tree with regular tree object.");
            }
            
            SeedbankTreeX other = (SeedbankTreeX)m_initial.get();
            root = other.root.copy();
            nodeCount = other.nodeCount;
            internalNodeCount = other.internalNodeCount;
            leafNodeCount = other.leafNodeCount;
        }
		// TODO: not sure why this is done
        if (nodeCount < 0) {
            if (m_taxonset.get() != null) {
                // make a caterpillar
                List<String> sTaxa = m_taxonset.get().asStringList();
                Node left = new SeedbankNodeX();
                left.setNr(0);
                left.setHeight(0);
                left.setID(sTaxa.get(0));
                for (int i = 1; i < sTaxa.size(); i++) {
                    Node right = new SeedbankNodeX();
                    right.setNr(i);
                    right.setHeight(0);
                    right.setID(sTaxa.get(i));
                    Node parent = new SeedbankNodeX();
                    parent.setNr(sTaxa.size() + i - 1);
                    parent.setHeight(i);
                    left.setParent(parent);
                    parent.setLeft(left);
                    right.setParent(parent);
                    parent.setRight(right);
                    left = parent;
                }
                root = left;
                leafNodeCount = sTaxa.size();
                nodeCount = leafNodeCount * 2 - 1;
                internalNodeCount = leafNodeCount - 1;

            } else {
                // make dummy tree with a single root node
                root = new SeedbankNodeX();
                root.setNr(0);
                root.setTree(this);
                nodeCount = 1;
                internalNodeCount = 0;
                leafNodeCount = 1;
            }
        }

        if (nodeCount >= 0) {
            initArrays();
        }
        
        typeLabel = typeLabelInput.get();
        
        processTraits(m_traitList.get());

        // Ensure tree is compatible with traits.
        if (hasDateTrait())
            adjustTreeNodeHeights(root);
        
     
     // SCMTT
     // Obtain required parameters from inputs:
        transitionModel = transitionModelInput.get();
        
     // Obtain leaf colours from explicit input or alignment:
        leafTypes = Lists.newArrayList();
        leafNames = Lists.newArrayList();
        if (leafTypesInput.get() != null) {
            for (int i=0; i<leafTypesInput.get().getDimension(); i++) {
                leafTypes.add(leafTypesInput.get().getValue(i));
                leafNames.add(String.valueOf(i));
            }
        } else {
            // Fill leaf colour array:
            if (hasTypeTrait()) {
                for (int i = 0; i<typeTraitSet.taxaInput.get().asStringList().size(); i++) {
//                	leafTypes.add(typeTraitSet.getValue(typeTraitSet.taxaInput.get().asStringList().get(i));
                    leafTypes.add(migModel.getTypeSet().getTypeIndex((typeTraitSet.getStringValue(i))));
                    leafNames.add(typeTraitSet.taxaInput.get().asStringList().get(i));
                }
            } else {
                throw new IllegalArgumentException("Either leafColours or "
                    + "trait set (with name '" + typeLabel + "') "
                    + "must be provided.");
            }

        }

//        // Count unique leaf types:
//        int nUniqueLeafTypes = new HashSet<>(leafTypes).size();
//        if (nUniqueLeafTypes > migModel.getNTypes())
//            throw new IllegalArgumentException("There are " + nUniqueLeafTypes
//                    + " unique leaf types but the model only includes "
//                    + migModel.getNTypes() + " unique types!");

        nLeaves = leafTypes.size();
        
        // Set leaf times if specified:
        leafTimes = Lists.newArrayList();
        if (timeTraitSet == null) {
            for (int i=0; i<nLeaves; i++)
                leafTimes.add(0.0);
        } else {
            if (timeTraitSet.taxaInput.get().asStringList().size() != nLeaves)
                throw new IllegalArgumentException("Number of time traits "
                        + "doesn't match number of leaf colours supplied.");
            
            for (int i=0; i<nLeaves; i++)
                leafTimes.add(timeTraitSet.getValue(i));
        }
        

        // Construct tree
        this.root = simulateTree();
        this.root.setParent(null);
        this.nodeCount = this.root.getNodeCount();
        this.internalNodeCount = this.root.getInternalNodeCount();
        this.leafNodeCount = this.root.getLeafNodeCount();
        initArrays();

        // Write tree to disk if requested
        if (outputFileNameInput.get() != null) {
            try (PrintStream pstream = new PrintStream(outputFileNameInput.get())) {
                pstream.println("#nexus\nbegin trees;");
                pstream.println("tree TREE_1 = " + toString() + ";");
                pstream.println("end;");
            } catch (FileNotFoundException e) {
                throw new RuntimeException("Error opening file '"
                        + outputFileNameInput.get() + "' for writing.");
            }
        }
	}
	
    @Override
    protected void processTraits(List<TraitSet> traitList) {
        super.processTraits(traitList);
        
        // Record trait set associated with leaf types.
        for (TraitSet traitSet : traitList) {
            if (traitSet.getTraitName().equals(typeLabel)) {
                typeTraitSet = traitSet;
                break;
            }
        }

        // Use explicitly-identified type trait set if available.
        // Seems dumb, but needed for BEAUti as ListInputEditors
        // muck things up...
        if (typeTraitInput.get() != null)
            typeTraitSet = typeTraitInput.get();

        // Construct type list.
        if (typeTraitSet == null) {
            if (getTaxonset() != null) {
                TraitSet dummyTraitSet = new TraitSet();

                StringBuilder sb = new StringBuilder();
                for (int i=0; i<getTaxonset().getTaxonCount(); i++) {
                    if (i>0)
                        sb.append(",\n");
                    sb.append(getTaxonset().getTaxonId(i)).append("=NOT_SET");
                }
                try {
                    dummyTraitSet.initByName(
                        "traitname", "type",
                        "taxa", getTaxonset(),
                        "value", sb.toString());
                    dummyTraitSet.setID("typeTraitSetInput.t:"
                        + BeautiDoc.parsePartition(getID()));
                    setTypeTrait(dummyTraitSet);

//                    if (typeSet != null)
//                        typeSet.addTypesFromTypeTraitSet(dummyTraitSet);
                } catch (Exception ex) {
                    System.out.println("Error setting default type trait.");
                }
            }
        }

//        if (typeSet == null) {
//            typeSet = new TypeSet();
//            typeSet.initByName("typeTraitSet", typeTraitSet);
//        }
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
    
    @Override
    public final void initArrays() {
        // initialise tree-as-array representation + its stored variant
        m_nodes = new SeedbankNodeX[nodeCount];
        listNodes((SeedbankNodeX)root, (SeedbankNodeX[])m_nodes);
        m_storedNodes = new SeedbankNodeX[nodeCount];
        Node copy = root.copy();
        listNodes((SeedbankNodeX)copy, (SeedbankNodeX[])m_storedNodes);
    }
    
    /**
     * Convert seedbank tree to array representation.
     *
     * @param node Root of sub-tree to convert.
     * @param nodes Array to populate with tree nodes.
     */
    private void listNodes(SeedbankNodeX node, SeedbankNodeX[] nodes) {
        nodes[node.getNr()] = node;
        node.setTree(this);
        if (!node.isLeaf()) {
            listNodes((SeedbankNodeX)node.getLeft(), nodes);
            if (node.getRight()!=null)
                listNodes((SeedbankNodeX)node.getRight(), nodes);
        }
    }
    
    /**
     *      * Deep copy, returns a completely new seedbank tree.
     *
     * @return a deep copy of this seedbank tree
     */
    @Override
    public SeedbankTreeX copy() {
    	SeedbankTreeX tree = new SeedbankTreeX();
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
        SeedbankTreeX mtTree = (SeedbankTreeX) other;

        SeedbankNodeX[] mtNodes = new SeedbankNodeX[mtTree.getNodeCount()];
        for (int i=0; i<mtTree.getNodeCount(); i++)
            mtNodes[i] = new SeedbankNodeX();

        ID = mtTree.ID;
        root = mtNodes[mtTree.root.getNr()];
        root.assignFrom(mtNodes, mtTree.root);
        root.setParent(null);

        nodeCount = mtTree.nodeCount;
        internalNodeCount = mtTree.internalNodeCount;
        leafNodeCount = mtTree.leafNodeCount;
        initArrays();
    }
    
    /**
     * Copy all values aside from IDs from an existing seedbank tree.
     * 
     * @param other
     */
    @Override
    public void assignFromFragile(StateNode other) {
    	SeedbankTreeX mtTree = (SeedbankTreeX) other;
        if (m_nodes == null) {
            initArrays();
        }
        root = m_nodes[mtTree.root.getNr()];
        Node[] otherNodes = mtTree.m_nodes;
        int iRoot = root.getNr();
        assignFromFragileHelper(0, iRoot, otherNodes);
        
        root.setHeight(otherNodes[iRoot].getHeight());
        root.setParent(null);
        
        SeedbankNodeX mtRoot = (SeedbankNodeX)root;
        mtRoot.nodeType = ((SeedbankNodeX)(otherNodes[iRoot])).nodeType;
//        mtRoot.changeTimes.clear();
//        mtRoot.changeTypes.clear();
//        mtRoot.nTypeChanges = 0;
        
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
     * helper to assignFromFragile *
     */
    private void assignFromFragileHelper(int iStart, int iEnd, Node[] otherNodes) {
        for (int i = iStart; i < iEnd; i++) {
        	SeedbankNodeX sink = (SeedbankNodeX)m_nodes[i];
        	SeedbankNodeX src = (SeedbankNodeX)otherNodes[i];
            
            sink.setHeight(src.getHeight());
            sink.setParent(m_nodes[src.getParent().getNr()]);
            
//            sink.nTypeChanges = src.nTypeChanges;
//            sink.changeTimes.clear();
//            sink.changeTimes.addAll(src.changeTimes);
//            sink.changeTypes.clear();
//            sink.changeTypes.addAll(src.changeTypes);
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
    
    // check that changes decrease in time from parent to child
    private boolean timesAreValid(Node node) {
        for (Node child : node.getChildren()) {
            if (child.getHeight()>node.getHeight())
                return false;
            
            if (!timesAreValid(child))
                return false;
        }
        
        return true;
    }

    // TODO: should there be a check that a type change is actually a change? 
    private boolean typesAreValid(Node node) {
        for (Node child : node.getChildren()) {
        	if (((SeedbankNodeX)node).getNodeEvent() == NodeEvent.MIGRATION) {
        		if (((SeedbankNodeX)node).getNodeToType() != ((SeedbankNodeX)child).getNodeType()) {
        			return false;
        		}
        	} else {
        		if (((SeedbankNodeX)node).getNodeType() != ((SeedbankNodeX)child).getNodeType()) {
        			return false;
        		}
        	}
            
            if (!typesAreValid(child))
                return false;
        }
        return true;
    }
    
    // SCMTT Methods
    
    /**
     * Generates tree using the specified list of active leaf nodes using the
     * structured coalescent.
     *
     * @return Root node of generated tree.
     */
    private SeedbankNodeX simulateTree() {

        // Initialise node creation counter:
        int nextNodeNr = 0;

        // Initialise node lists:
        List<List<SeedbankNodeX>> liveNodes = Lists.newArrayList();
        List<List<SeedbankNodeX>> deadNodes = Lists.newArrayList();

        liveNodes.add(new ArrayList<>());
        liveNodes.add(new ArrayList<>());
        deadNodes.add(new ArrayList<>());
        deadNodes.add(new ArrayList<>());

        // Add leaves to dead nodes list:
        for (int l = 0; l < nLeaves; l++) {
            SeedbankNodeX node = new SeedbankNodeX();
            node.setNr(nextNodeNr);
            node.setID(leafNames.get(l));
            deadNodes.get(leafTypes.get(l)).add(node);
            node.setHeight(leafTimes.get(l));
            node.setNodeType(leafTypes.get(l));
            node.setNodeEvent(NodeEvent.SAMPLE);

            nextNodeNr++;
        }
        
        // Sort nodes in dead nodes lists in order of increasing age:
        for (int i=0; i<2; i++) {
            Collections.sort(deadNodes.get(i),
                (SeedbankNodeX node1, SeedbankNodeX node2) -> {
                    double dt = node1.getHeight()-node2.getHeight();
                    if (dt<0)
                        return -1;
                    if (dt>0)
                        return 1;
                    
                    return 0;
                });
        }

        // Allocate propensity lists:
        List<Double> migrationProp = new ArrayList<>();
        Double coalesceProp = 0.0;
        double t = 0.0;

        while (totalNodesRemaining(liveNodes)>1
                || totalNodesRemaining(deadNodes)>0) {

            // Step 1: Calculate propensities.
            double totalProp = updatePropensities(migrationProp, coalesceProp,
                    liveNodes);
            
            // Step 2: Determine next event.
            SBEvent event = getNextEvent(migrationProp, coalesceProp,
                    totalProp, t);

            // Step 3: Handle activation of nodes:
            SeedbankNodeX nextNode = null;
            int nextNodeType = -1;
            double nextTime = Double.POSITIVE_INFINITY;
            
            for (int i=0; i<2; i++) {
                if (deadNodes.get(i).isEmpty())
                    continue;
                
                if (deadNodes.get(i).get(0).getHeight()<nextTime) { // finds the youngest deadNode
                    nextNode = deadNodes.get(i).get(0);
                    nextTime = nextNode.getHeight();
                    nextNodeType = i;
                }
            }
            if (nextTime < event.time) { // if birth time sooner than event time
                t = nextTime;
                liveNodes.get(nextNodeType).add(nextNode);
                deadNodes.get(nextNodeType).remove(0);
                continue;
            }
            
            // Step 4: Place event on tree.
            nextNodeNr = updateTree(liveNodes, event, nextNodeNr);

            // Step 5: Keep track of time increment.
            t = event.time;
        }

        // TODO: assert here that the remaining live node must be of the active type?
        // Return sole remaining live node as root:
        for (List<SeedbankNodeX> nodeList : liveNodes)
            if (!nodeList.isEmpty())
                return nodeList.get(0);

        // Should not fall through.
        throw new RuntimeException("No live nodes remaining end of "
                + "structured coalescent simulation!");
    }
    
    /**
     * Obtain propensities (instantaneous reaction rates) for coalescence and
     * migration events.
     *
     * @param migrationProp
     * @param coalesceProp
     * @param liveNodes
     * @return Total reaction propensity.
     */
    private double updatePropensities(List<Double> migrationProp,
            Double coalesceProp, List<List<SeedbankNodeX>> liveNodes) {

        double totalProp = 0.0;
        
        double N_a = transitionModel.getPopSize(1);
        int k_a = liveNodes.get(1).size();
        int k_d = liveNodes.get(0).size();
        double m_ad = transitionModel.getBackwardRate(1, 0);
        double m_da = transitionModel.getBackwardRate(0, 1);
        
        coalesceProp = k_a * (k_a - 1) / (2.0 * N_a);
        totalProp += coalesceProp;
        
        migrationProp.set(1, k_a * m_ad);
        migrationProp.set(0, k_d * m_da);
        totalProp += migrationProp.get(1);
        totalProp += migrationProp.get(0);

        return totalProp;

    }
    
    /**
     * Calculate total number of live nodes remaining.
     *
     * @param liveNodes
     * @return Number of live nodes remaining.
     */
    private int totalNodesRemaining(List<List<SeedbankNodeX>> liveNodes) {
        return liveNodes.get(0).size() + liveNodes.get(1).size();
    }
    
    /**
     * Obtain type and location of next reaction.
     *
     * @param migrateProp Current migration propensities.
     * @param coalesceProp Current coalescence propensities.
     * @param t Current time.
     * @return Event object describing next event.
     */
    private SBEvent getNextEvent(List<Double> migrateProp,
            Double coalesceProp, double totalProp, double t) {

        // Get time of next event:
        if (totalProp>0.0)
            t += Randomizer.nextExponential(totalProp);
        else
            return new NullEvent();

        // Select event type:
        double U = Randomizer.nextDouble() * totalProp;
        
        if (U < coalesceProp) //TODO
        	return new CoalescenceEvent(1, t);
        else
        	U -= coalesceProp;
        
        if (U < migrateProp.get(0))
            return new MigrationEvent(0, 1, t);
        else
            U -= migrateProp.get(0);
        
        if (U < migrateProp.get(1))
            return new MigrationEvent(1, 0, t);
        else
            U -= migrateProp.get(1);

        // Should not fall through.
        throw new RuntimeException("Structured coalescenct event selection error.");
    }
    
    /**
    * Update tree with result of latest event.
    *
    * @param liveNodes
    * @param event
    * @param nextNodeNr Integer identifier of last node added to tree.
    * @return Updated nextNodeNr.
    */
   private int updateTree(List<List<SeedbankNodeX>> liveNodes, SBEvent event,
           int nextNodeNr) {

       if (event instanceof CoalescenceEvent) {

           // Randomly select node pair from active nodes:
    	   SeedbankNodeX daughter = selectRandomNode(liveNodes.get(1));
    	   SeedbankNodeX son = selectRandomSibling(liveNodes.get(1), daughter);

           // Create new parent node with appropriate ID and time:
    	   SeedbankNodeX parent = new SeedbankNodeX();
           parent.setNr(nextNodeNr);
           parent.setID(String.valueOf(nextNodeNr));
           parent.setHeight(event.time);
           parent.setNodeEvent(NodeEvent.COALESCENT);
           nextNodeNr++;

           // Connect new parent to children:
           parent.setLeft(daughter);
           parent.setRight(son);
           son.setParent(parent);
           daughter.setParent(parent);

           // Ensure new parent is set to correct colour:
           parent.setNodeType(event.fromType);

           // Update activeNodes:
           liveNodes.get(1).remove(son);
           int idx = liveNodes.get(1).indexOf(daughter);
           liveNodes.get(1).set(idx, parent); // update daughter idx to parent

       } else { // Migration event ... what happens if Null event?

           // Randomly select node with chosen colour:
    	   SeedbankNodeX daughter = selectRandomNode(liveNodes.get(event.fromType));
    	   
    	   // // !!
    	   SeedbankNodeX migrator = new SeedbankNodeX();

           migrator.setNr(nextNodeNr);
           migrator.setID(String.valueOf(nextNodeNr));
           migrator.setHeight(event.time);
           migrator.setNodeEvent(NodeEvent.MIGRATION);
           nextNodeNr++;
           
           // SeedbankNode ToType is time forward
           migrator.setNodeType(event.toType);
           migrator.setNodeToType(event.fromType);
           
           // Make dummy right child
           SeedbankNodeX dummy = new SeedbankNodeX();
           dummy.setNr(nextNodeNr);
           dummy.setID(String.valueOf(nextNodeNr));
           dummy.setHeight(event.time); // branch length 0
           dummy.setNodeEvent(NodeEvent.DUMMY);
           nextNodeNr++;
           
           dummy.setNodeType(event.fromType);

           // Connect new parent to children:
           migrator.setLeft(daughter);
           migrator.setRight(dummy);
           dummy.setParent(migrator);
           daughter.setParent(migrator);


           // Update activeNodes:
           liveNodes.get(event.fromType).remove(daughter);
           liveNodes.get(event.toType).add(migrator);

       }

       return nextNodeNr;

   }
    
    /**
     * Use beast RNG to select random node from list.
     *
     * @param nodeList
     * @return A randomly selected node.
     */
    private SeedbankNodeX selectRandomNode(List<SeedbankNodeX> nodeList) {
        return nodeList.get(Randomizer.nextInt(nodeList.size()));
    }

    /**
     * Return random node from list, excluding given node.
     *
     * @param nodeList
     * @param node
     * @return Randomly selected node.
     */
    private SeedbankNodeX selectRandomSibling(List<SeedbankNodeX> nodeList, Node node) {

        int n = Randomizer.nextInt(nodeList.size() - 1);
        int idxToAvoid = nodeList.indexOf(node);
        if (n >= idxToAvoid)
            n++;

        return nodeList.get(n);
    }
    
    @Override
    public void initStateNodes() { }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodeList) {
        stateNodeList.add(this);
    }
    
    /**
     * Generates an ensemble of trees from the structured coalescent for testing
     * coloured tree-space samplers.
     *
     * @param argv
     * @throws java.lang.Exception
     */
    public static void main(String[] argv) throws Exception {
        
//        // Set up migration model.
//        RealParameter rateMatrix = new RealParameter();
//        rateMatrix.initByName(
//                "value", "0.05",
//                "dimension", "12");
//        RealParameter popSizes = new RealParameter();
//        popSizes.initByName(
//                "value", "7.0",
//                "dimension", "4");
//        SCMigrationModel migrationModel = new SCMigrationModel();
//        migrationModel.initByName(
//                "rateMatrix", rateMatrix,
//                "popSizes", popSizes);
//
//        // Specify leaf types:
//        IntegerParameter leafTypes = new IntegerParameter();
//        leafTypes.initByName(
//                "value", "0 0 0");
//
//        // Generate ensemble:
//        int reps = 1000000;
//        double[] heights = new double[reps];
//        double[] changes = new double[reps];
//
//        long startTime = System.currentTimeMillis();
//        StructuredCoalescentMultiTypeTree sctree;
//        sctree = new StructuredCoalescentMultiTypeTree();
//        for (int i = 0; i < reps; i++) {
//
//            if (i % 1000 == 0)
//                System.out.format("%d reps done\n", i);
//
//            sctree.initByName(
//                    "migrationModel", migrationModel,
//                    "leafTypes", leafTypes,
//                    "nTypes", 4);
//
//            heights[i] = sctree.getRoot().getHeight();
//            changes[i] = sctree.getTotalNumberOfChanges();
//        }
//
//        long time = System.currentTimeMillis() - startTime;
//
//        System.out.printf("E[T] = %1.4f +/- %1.4f\n",
//                DiscreteStatistics.mean(heights), DiscreteStatistics.stdev(heights) / Math.sqrt(reps));
//        System.out.printf("V[T] = %1.4f\n", DiscreteStatistics.variance(heights));
//
//        System.out.printf("E[C] = %1.4f +/- %1.4f\n",
//                DiscreteStatistics.mean(changes), DiscreteStatistics.stdev(changes) / Math.sqrt(reps));
//        System.out.printf("V[C] = %1.4f\n", DiscreteStatistics.variance(changes));
//
//        System.out.printf("Took %1.2f seconds\n", time / 1000.0);
//
//        try (PrintStream outStream = new PrintStream("heights.txt")) {
//            outStream.println("h c");
//            for (int i = 0; i < reps; i++)
//                outStream.format("%g %g\n", heights[i], changes[i]);
//        }
    }
    
}
