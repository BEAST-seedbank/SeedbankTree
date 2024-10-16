 package seedbanktree.distributions;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import seedbanktree.evolution.tree.SeedbankNode;
import seedbanktree.evolution.tree.SeedbankTree;
import seedbanktree.evolution.tree.TransitionModel;

public class SeedbankTreeDensity extends Distribution {
	
	public Input<SeedbankTree> sbTreeInput = new Input<>("tree",
			"Seedbank tree.", Validate.REQUIRED);

	public Input<TransitionModel> transitionModelInput = new Input<>(
            "transitionModel", "Model of transition between activity and dormancy.",
            Validate.REQUIRED);
	
    public Input<Boolean> checkValidityInput = new Input<>(
            "checkValidity", "Explicitly check validity of colouring.  "
            +"(Default false.)  "
            +"Useful if operators are in danger of proposing invalid trees.",
            false);
    
    public Input<Boolean> debugLoggingInput = new Input<> (
    		"debugLogging", "Print logging statements for testing. (Default false.)", false);
    
    private SeedbankTree sbTree;
    private TransitionModel transitionModel;
    private boolean checkValidity, debugLogging;
    
    private enum SBEventKind {
        COALESCE, MIGRATE, SAMPLE
    };

    private class SBEvent {
        double time;
        int type, destType;
        SBEventKind kind;
        Node node;
    }
    
    private List<SBEvent> eventList;
    private List<Integer[]> lineageCountList;

	@Override
    public void initAndValidate() {
        sbTree = sbTreeInput.get();
        transitionModel = transitionModelInput.get();
        checkValidity = checkValidityInput.get();
        debugLogging = debugLoggingInput.get();
        
        eventList = new ArrayList<>();
        lineageCountList = new ArrayList<>();
	}

	@Override
    public double calculateLogP() {
        
        // Check validity of tree if required:
        if (checkValidity && sbTree.somethingIsDirty() && !sbTree.isValid())
            return Double.NEGATIVE_INFINITY;

        // Ensure sequence of events is up-to-date:
        updateEventSequence();

        // Start from the tips of the tree, working up.
        logP = 0;

        // Note that the first event is always a sample. We begin at the first
        // _interval_ and the event following that interval.
        for (int eventIdx = 1; eventIdx<eventList.size(); eventIdx++) {

            SBEvent event = eventList.get(eventIdx);
            Integer[] lineageCount = lineageCountList.get(eventIdx);
            double delta_t = event.time-eventList.get(eventIdx-1).time;

            // Interval contribution:
            if (delta_t>0) {
                double lambda = 0.0;
                double N_a = transitionModel.getPopSize(1);
                double theta = N_a * 2 * (1);
                theta = N_a;
                theta = 1;
                int k_a = lineageCount[1];
                int k_d = lineageCount[0];
                double m_ad = transitionModel.getBackwardRate(1, 0);
                double m_da = transitionModel.getBackwardRate(0, 1);
                
                
                lambda += k_a*(k_a-1)/(2.0*theta);
                lambda += k_a*m_ad;
                lambda += k_d*m_da;
                		
                logP += -delta_t*lambda;
                
                if (debugLogging) {
                	System.out.println("TIME CONTRIBUTION");
                	System.out.println(String.format("Interval: %f\nlogP: %f\n", delta_t, -delta_t*lambda));
                	manualLog("spikeandslab/manualLogging.txt", "TIME CONTRIBUTION");
                	manualLog("spikeandslab/manualLogging.txt", String.format("Interval: %f\nlogP: %f\n", delta_t, -delta_t*lambda));
                }
            }

            // Event contribution:
            switch (event.kind) {
                case COALESCE:
                    double N = transitionModel.getPopSize(event.type);
                    double theta = N * 2 * (1);
                    theta = N;
                    theta = 1;
                    logP += Math.log(1.0/theta);
                    
                    if (debugLogging) {
                    	System.out.print("COALESCE EVENT: ");
                        System.out.println(String.format("logP: %f", Math.log(1.0/theta)));
                        manualLog("spikeandslab/manualLogging.txt", "COALESCE EVENT: ");
                        manualLog("spikeandslab/manualLogging.txt", String.format("logP: %f", Math.log(1.0/theta)));
                    }

                    break;

                case MIGRATE:
                    double m = transitionModel
                            .getBackwardRate(event.type, event.destType);
                    logP += Math.log(m);
                    
                    if (debugLogging) {
                    	System.out.print(String.format("MIGRATE EVENT: %d to %d", event.type, event.destType));
                        System.out.println(String.format("logP: %f", Math.log(m)));
                        manualLog("spikeandslab/manualLogging.txt", String.format("MIGRATE EVENT: %d to %d", event.type, event.destType));
                        manualLog("spikeandslab/manualLogging.txt", String.format("logP: %f", Math.log(m)));
                    }

                    break;

                case SAMPLE:
                    // Do nothing here: only effect of sampling event is
                    // to change the lineage counts in subsequent intervals.
                	
                	if (debugLogging) {
                		System.out.print("SAMPLE: ");
                        System.out.println("logP: 0");
                        manualLog("spikeandslab/manualLogging.txt", "SAMPLE: ");
                        manualLog("spikeandslab/manualLogging.txt", "logP: 0");
                        
                	}
                	
                    break;
            }
            
            if (debugLogging) {
        		System.out.println("\nTotal logP: " + logP);
        		System.out.println("---");
        		manualLog("spikeandslab/manualLogging.txt", "\nTotal logP: " + logP);
        		manualLog("spikeandslab/manualLogging.txt", "---");
        		manualLog("spikeandslab/manualLogging.txt", "");
        		
        	}
        }
        
        return logP;
    }

    /**
     * Determines the sequence of migration, coalescence and sampling events
     * which make up the seedbank tree.
     */
    protected void updateEventSequence() {

        // Clean up previous list:
        eventList.clear();
        lineageCountList.clear();
        Node rootNode = sbTree.getRoot();

        // Initialise map of active nodes to active change indices:
        Map<Node, Integer> changeIdx = new HashMap<>();
        changeIdx.put(rootNode, -1);

        // Initialise lineage count:
        Integer[] lineageCount = new Integer[2];
        lineageCount[0] = ((SeedbankNode)rootNode).getNodeType() == 0 ? 1 : 0;
        lineageCount[1] = ((SeedbankNode)rootNode).getNodeType() == 1 ? 1 : 0;

        // Calculate event sequence:
        while (!changeIdx.isEmpty()) {

            SBEvent nextEvent = new SBEvent();
            nextEvent.time = Double.NEGATIVE_INFINITY;
            nextEvent.node = rootNode; // Initial assignment not significant

            // Determine next event
            for (Node node : changeIdx.keySet())
                if (changeIdx.get(node)<0) {
                    if (node.isLeaf()) {
                        // Next event is a sample
                        if (node.getHeight()>nextEvent.time) {
                            nextEvent.time = node.getHeight();
                            nextEvent.kind = SBEventKind.SAMPLE;
                            nextEvent.type = ((SeedbankNode)node).getNodeType();
                            nextEvent.node = node;
                        }
                    } else {
                        // Next event is a coalescence
                        if (node.getHeight()>nextEvent.time) {
                            nextEvent.time = node.getHeight();
                            nextEvent.kind = SBEventKind.COALESCE;
                            nextEvent.type = ((SeedbankNode)node).getNodeType();
                            nextEvent.node = node;
                        }
                    }
                } else {
                    // Next event is a migration
                    double thisChangeTime = ((SeedbankNode)node).getChangeTime(changeIdx.get(node));
                    if (thisChangeTime>nextEvent.time) {
                        nextEvent.time = thisChangeTime;
                        nextEvent.kind = SBEventKind.MIGRATE;
                        nextEvent.destType = ((SeedbankNode)node).getChangeType(changeIdx.get(node));
                        if (changeIdx.get(node)>0)
                            nextEvent.type = ((SeedbankNode)node).getChangeType(changeIdx.get(node)-1);
                        else
                            nextEvent.type = ((SeedbankNode)node).getNodeType();
                        nextEvent.node = node;
                    }
                }

            // Update active node list (changeIdx) and lineage count appropriately:
            switch (nextEvent.kind) {
                case COALESCE:
                    Node leftChild = nextEvent.node.getLeft();
                    Node rightChild = nextEvent.node.getRight();

                    changeIdx.remove(nextEvent.node);
                    changeIdx.put(leftChild, ((SeedbankNode)leftChild).getChangeCount()-1);
                    changeIdx.put(rightChild, ((SeedbankNode)rightChild).getChangeCount()-1);
                    lineageCount[nextEvent.type]++;
                    break;

                case SAMPLE:
                    changeIdx.remove(nextEvent.node);
                    lineageCount[nextEvent.type]--;
                    break;

                case MIGRATE:
                    lineageCount[nextEvent.destType]--;
                    lineageCount[nextEvent.type]++;
                    int oldIdx = changeIdx.get(nextEvent.node);
                    changeIdx.put(nextEvent.node, oldIdx-1);
                    break;
            }

            // Add event to list:
            eventList.add(nextEvent);
            lineageCountList.add(Arrays.copyOf(lineageCount, lineageCount.length));
        }

        // Reverse event and lineage count lists (order them from tips to root):
        Collections.reverse(eventList);
        Collections.reverse(lineageCountList);

    }

    @Override
    public boolean requiresRecalculation() {
    	// Assumption that if any inputs change, there will be recalculation
        return true;
    }
	
	
	// Distribution interface requirements
	
    /**
     * @return a list of unique ids for the state nodes that make up the arguments
     */
	@Override
	public List<String> getArguments() {
		return Collections.singletonList(sbTreeInput.get().getID());
	}

	/**
     * @return a list of unique ids for the state nodes that make up the conditions
     */
	@Override
	public List<String> getConditions() {
		return transitionModelInput.get().getParameterIds();
	}

	@Override
	public void sample(State state, Random random) {	
	}

}
