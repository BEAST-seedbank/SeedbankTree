package seedbanktree.evolution.tree;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;

public class SeedbankTreeInitialiser extends SeedbankTree implements StateNodeInitialiser {
	public Input<RealParameter> lambdasInput = new Input<> ("lambdas", "");
    public Input<IntegerParameter> indicatorsInput = new Input<> ("indicators", "");
    
	public Input<TransitionModel> transitionModelInput = new Input<>("transitionModel", "transition model to use in simulator.", Validate.REQUIRED);
	
	private TransitionModel transitionModel;
    
    private List<Integer> leafTypes;
    private List<String> leafNames;
    private List<Double> leafTimes;
    private int nLeaves;
    
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
    
    public SeedbankTreeInitialiser() {
    	m_initial.setRule(Validate.REQUIRED);
    }
    
    @Override
    public void initAndValidate() {
    	
    	super.initAndValidate();
    	
        // Obtain required parameters from inputs:
           transitionModel = transitionModelInput.get();
           
        // Obtain leaf colours from explicit input or alignment:
           leafTypes = new ArrayList<>();
           leafNames = new ArrayList<>();
        
   	    // Fill leaf colour array:
           assert(hasTypeTrait());
           for (int i = 0; i<typeTraitSet.taxaInput.get().asStringList().size(); i++) {
               leafTypes.add(getTypeIndex((typeTraitSet.getStringValue(i))));
               leafNames.add(typeTraitSet.taxaInput.get().asStringList().get(i));
           }

           nLeaves = leafTypes.size();
           
           // Set leaf times if specified:
           leafTimes = new ArrayList<>();
           if (!hasDateTrait()) {
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
           
           initStateNodes();
    }
    
    /**
     * Generates tree using the specified list of active leaf nodes using the
     * seedbank coalescent.
     *
     * @return Root node of generated tree.
     */
    private SeedbankNode simulateTree() {

        // Initialise node creation counter:
        int nextNodeNr = 0;

        // Initialise node lists:
        List<List<SeedbankNode>> liveNodes = new ArrayList<>();
        List<List<SeedbankNode>> deadNodes = new ArrayList<>();
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
        migrationProp.add(0.0);
        migrationProp.add(0.0);
        List<Double> coalesceProp = new ArrayList<>();
        coalesceProp.add(0.0);
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
        
        // Return sole remaining live node as root:
        for (List<SeedbankNode> nodeList : liveNodes)
            if (!nodeList.isEmpty()) {
            	((SeedbankNode) nodeList.get(0).getLeft()).clearChanges();
            	return nodeList.get(0);
            }

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
            List<Double> coalesceProp, List<List<SeedbankNode>> liveNodes) {

        double totalProp = 0.0;
        
        double N_a = transitionModel.getPopSize(1);
        int k_a = liveNodes.get(1).size();
        int k_d = liveNodes.get(0).size();
        double m_ad = transitionModel.getBackwardRate(1, 0);
        double m_da = transitionModel.getBackwardRate(0, 1);
        
        coalesceProp.set(0, k_a * (k_a - 1) / (2.0 * N_a));
        coalesceProp.set(0, k_a * (k_a - 1) / (2.0 * 1));
        totalProp += coalesceProp.get(0);
        
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
            List<Double> coalesceProp, double totalProp, double t) {
        // Get time of next event:
        if (totalProp>0.0)
            t += Randomizer.nextExponential(totalProp);
        else
            return new NullEvent();

        // Select event type:
        double U = Randomizer.nextDouble() * totalProp;
        
        if (U < coalesceProp.get(0)) 
        	return new CoalescenceEvent(1, t);
        else
        	U -= coalesceProp.get(0);
        
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
        assert(m_initial != null);
        if (!(m_initial.get() instanceof SeedbankTree)) {
        	throw new IllegalArgumentException("Attempted to use "
                    + "seedbank tree initialiser on regular tree object");
        }
    	m_initial.get().assignFromWithoutID(this);
    	
    	if (lambdasInput.get() == null || indicatorsInput.get() == null) {
    		return;
    	}
    				
    	RealParameter lambdas = lambdasInput.get();
        IntegerParameter indicators = indicatorsInput.get();
        for (int i = 0; i < getNodeCount()-1; i++) {
        	SeedbankNode sbNode = (SeedbankNode)getNode(i);
        	if (sbNode.getChangeCount() != 0) {
        		double sum = 0.0;
        		int lastType = sbNode.getNodeType();
        		double lastHeight = sbNode.getHeight();
        		for (int j = 0; j < sbNode.getChangeCount(); j++) {
        			if (lastType == 0)
        				sum += sbNode.getChangeTime(j) - lastHeight;
        			lastType = sbNode.getChangeType(j);
        			lastHeight = sbNode.getChangeTime(j);
        		}
        		indicators.setValue(i, 1);
        		lambdas.setValue(i, sum / sbNode.getLength());
        	}
        }
    }

    @Override
    public void getInitialisedStateNodes(final List<StateNode> stateNodes) {
    	assert(m_initial != null);
        if (!(m_initial.get() instanceof SeedbankTree)) {
        	throw new IllegalArgumentException("Attempted to use "
                    + "seedbank tree initialiser on regular tree object");
        }
        stateNodes.add(m_initial.get());
        stateNodes.add(lambdasInput.get());
        stateNodes.add(indicatorsInput.get());
    }

}
