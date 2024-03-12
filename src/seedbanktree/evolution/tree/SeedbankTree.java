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


public class SeedbankTree extends Tree implements StateNodeInitialiser {
	
	// Fields
	
	// Inputs MTT
	public Input<String> typeLabelInput = new Input<>("typeLabel", "Label for type traits (default 'type')", "type");
	
    public Input<TraitSet> typeTraitInput = new Input<>("typeTrait", "Type trait set.  Used only by BEAUti.");
    
    public Input<String> activeTypeNameInput = new Input<>(
            "activeTypeName", "Name of active type.", "active", Validate.OPTIONAL);
    
    public Input<String> dormantTypeNameInput = new Input<>(
            "activeTypeName", "Name of active type.", "active", Validate.OPTIONAL);
    
    // Inputs SCMTT
    public Input<TransitionModel> transitionModelInput = new Input<>(
            "transitionModel,",
            "Transition model to use in simulator.",
            Validate.REQUIRED);
    
	// Shadow inputs MTT
	private String typeLabel;
    private TraitSet typeTraitSet;
    private String activeTypeName;
    private String dormantTypeName;
    
    //SCMTT
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
            this.fromType = type;
            this.toType = type;
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
    
	
	@Override
	public void initAndValidate() {
		// MTT
		// If an initial tree is given as input
		if (m_initial.get() != null) {
            
            if (!(m_initial.get() instanceof SeedbankTree)) {
                throw new IllegalArgumentException("Attempted to initialise "
                        + "seedbank tree with regular tree object.");
            }
            
            SeedbankTree other = (SeedbankTree)m_initial.get();
            root = other.root.copy();
            nodeCount = other.nodeCount;
            internalNodeCount = other.internalNodeCount;
            leafNodeCount = other.leafNodeCount;
        }
        
        typeLabel = typeLabelInput.get();
        activeTypeName = activeTypeNameInput.get();
        dormantTypeName = dormantTypeNameInput.get();
        
        processTraits(m_traitList.get());
     
     // SCMTT
     // Obtain required parameters from inputs:
        transitionModel = transitionModelInput.get();
        
     // Obtain leaf colours from explicit input or alignment:
        leafTypes = Lists.newArrayList();
        leafNames = Lists.newArrayList();
     
	    // Fill leaf colour array:
        assert(hasTypeTrait());
        for (int i = 0; i<typeTraitSet.taxaInput.get().asStringList().size(); i++) {
            leafTypes.add(transitionModel.getTypeIndex((typeTraitSet.getStringValue(i))));
            leafNames.add(typeTraitSet.taxaInput.get().asStringList().get(i));
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
        SeedbankTree mtTree = (SeedbankTree) other;

        SeedbankNode[] mtNodes = new SeedbankNode[mtTree.getNodeCount()];
        for (int i=0; i<mtTree.getNodeCount(); i++)
            mtNodes[i] = new SeedbankNode();

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
    	SeedbankTree mtTree = (SeedbankTree) other;
        if (m_nodes == null) {
            initArrays();
        }
        root = m_nodes[mtTree.root.getNr()];
        Node[] otherNodes = mtTree.m_nodes;
        int iRoot = root.getNr();
        assignFromFragileHelper(0, iRoot, otherNodes);
        
        root.setHeight(otherNodes[iRoot].getHeight());
        root.setParent(null);
        
        SeedbankNode mtRoot = (SeedbankNode)root;
        mtRoot.nodeType = ((SeedbankNode)(otherNodes[iRoot])).nodeType;
        mtRoot.changeTimes.clear();
        mtRoot.changeTypes.clear();
        mtRoot.nTypeChanges = 0;
        
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
    
    // check that changes decrease in time from parent to child
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
        for (Node child : node.getChildren()) {
            if (((SeedbankNode)node).getNodeType() != ((SeedbankNode)child).getFinalType())
                return false;
            
            if (!typesAreValid(child))
                return false;
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
    
    // SCMTT Methods
    
    /**
     * Generates tree using the specified list of active leaf nodes using the
     * structured coalescent.
     *
     * @return Root node of generated tree.
     */
    private SeedbankNode simulateTree() {

        // Initialise node creation counter:
        int nextNodeNr = 0;

        // Initialise node lists:
        List<List<SeedbankNode>> liveNodes = Lists.newArrayList();
        List<List<SeedbankNode>> deadNodes = Lists.newArrayList();
//        for (int i = 0; i < migModel.getNTypes(); i++) {
//            activeNodes.add(new ArrayList<>());
//            inactiveNodes.add(new ArrayList<>());
//        }
        liveNodes.add(new ArrayList<>());
        liveNodes.add(new ArrayList<>());
        deadNodes.add(new ArrayList<>());
        deadNodes.add(new ArrayList<>());

        // Add nodes to dead nodes list:
        for (int l = 0; l < nLeaves; l++) {
            SeedbankNode node = new SeedbankNode();
            node.setNr(nextNodeNr);
            node.setID(leafNames.get(l));
            deadNodes.get(leafTypes.get(l)).add(node);
            node.setHeight(leafTimes.get(l));
            node.setNodeType(leafTypes.get(l));

            nextNodeNr++;
        }
        
        // Sort nodes in dead nodes lists in order of increasing age:
//        for (int i=0; i<migModel.getNTypes(); i++) {
        for (int i=0; i<2; i++) {
            Collections.sort(deadNodes.get(i),
                (SeedbankNode node1, SeedbankNode node2) -> {
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
            SeedbankNode nextNode = null;
            int nextNodeType = -1;
            double nextTime = Double.POSITIVE_INFINITY;
//            for (int i=0; i<migModel.getNTypes(); i++) {
            for (int i=0; i<2; i++) {
                if (deadNodes.get(i).isEmpty())
                    continue;
                
                if (deadNodes.get(i).get(0).getHeight()<nextTime) {
                    nextNode = deadNodes.get(i).get(0);
                    nextTime = nextNode.getHeight();
                    nextNodeType = i;
                }
            }
            if (nextTime < event.time) {
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
        for (List<SeedbankNode> nodeList : liveNodes)
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
     * @param activeNodes
     * @return Total reaction propensity.
     */
    private double updatePropensities(List<Double> migrationProp,
            Double coalesceProp, List<List<SeedbankNode>> liveNodes) {

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
    private int totalNodesRemaining(List<List<SeedbankNode>> liveNodes) {
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
        
        if (U < coalesceProp) 
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
   private int updateTree(List<List<SeedbankNode>> liveNodes, SBEvent event,
           int nextNodeNr) {

       if (event instanceof CoalescenceEvent) {

           // Randomly select node pair from active nodes:
    	   SeedbankNode daughter = selectRandomNode(liveNodes.get(1));
    	   SeedbankNode son = selectRandomSibling(liveNodes.get(1), daughter);

           // Create new parent node with appropriate ID and time:
    	   SeedbankNode parent = new SeedbankNode();
           parent.setNr(nextNodeNr);
           parent.setID(String.valueOf(nextNodeNr));
           parent.setHeight(event.time);
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
           liveNodes.get(1).set(idx, parent);

       } else { // Migration event ... what happens if Null event?

           // Randomly select node with chosen colour:
    	   SeedbankNode migrator = selectRandomNode(liveNodes.get(event.fromType));

           // Record colour change in change lists:
           migrator.addChange(event.toType, event.time);

           // Update activeNodes:
           liveNodes.get(event.fromType).remove(migrator);
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
    private SeedbankNode selectRandomNode(List<SeedbankNode> nodeList) {
        return nodeList.get(Randomizer.nextInt(nodeList.size()));
    }

    /**
     * Return random node from list, excluding given node.
     *
     * @param nodeList
     * @param node
     * @return Randomly selected node.
     */
    private SeedbankNode selectRandomSibling(List<SeedbankNode> nodeList, Node node) {

        int n = Randomizer.nextInt(nodeList.size() - 1);
        int idxToAvoid = nodeList.indexOf(node);
        if (n >= idxToAvoid)
            n++;

        return nodeList.get(n);
    }
    
    // Methods for StateNodeInitialiser interface
    
    @Override
    public void initStateNodes() {
        if (m_initial.get() != null) {
            m_initial.get().assignFromWithoutID(this);
        }
    }

    @Override
    public void getInitialisedStateNodes(final List<StateNode> stateNodes) {
        if (m_initial.get() != null) {
            stateNodes.add(m_initial.get());
        }
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
