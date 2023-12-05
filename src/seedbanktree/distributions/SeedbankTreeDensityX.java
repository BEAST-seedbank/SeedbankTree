package seedbanktree.distributions;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import com.google.common.collect.Lists;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import seedbanktree.evolution.tree.SeedbankNodeX;
import seedbanktree.evolution.tree.SeedbankNodeX.NodeEvent;
import seedbanktree.evolution.tree.SeedbankTreeX;
import seedbanktree.evolution.tree.TransitionModel;

public class SeedbankTreeDensityX extends Distribution{
	
	public Input<SeedbankTreeX> sbTreeInput = new Input<>("SeedbankTree",
			"Multi-type tree.", Validate.REQUIRED);

	public Input<TransitionModel> transitionModelInput = new Input<>(
            "transitionModel", "Model of transition between activity and dormancy.",
            Validate.REQUIRED);
	
    public Input<Boolean> checkValidityInput = new Input<>(
            "checkValidity", "Explicitly check validity of colouring.  "
            +"(Default false.)  "
            +"Useful if operators are in danger of proposing invalid trees.",
            false);
    
    protected SeedbankTreeX sbTree;
    protected TransitionModel transitionModel;
    protected boolean checkValidity;
    
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

    // Empty constructor as required:
    public SeedbankTreeDensityX() { };
    
    
	@Override
    public void initAndValidate() {
        sbTree = sbTreeInput.get();
        transitionModel = transitionModelInput.get();
        checkValidity = checkValidityInput.get();
        
        eventList = new ArrayList<>();
        lineageCountList = new ArrayList<>();
        
        // Ensure tree and migration model are compatible
//        if (mtTree.hasTypeTrait() && !mtTree.getTypeSet().equals(migrationModel.getTypeSet()))
//            throw new IllegalArgumentException("Tree and migration model have incompatible type sets.");
    }

	@Override
    public double calculateLogP() {
        
        // Check validity of tree if required:
        if (checkValidity && !sbTree.isValid())
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
                int k_a = lineageCount[1];
                int k_d = lineageCount[0];
                double m_ad = transitionModel.getBackwardRate(1, 0);
                double m_da = transitionModel.getBackwardRate(0, 1);
                
                
                lambda += k_a*(k_a-1)/(2.0*N_a);
                lambda += k_a*m_ad;
                lambda += k_d*m_da;
                		
                logP += -delta_t*lambda;
            }

            // Event contribution:
            switch (event.kind) {
                case COALESCE:
                    double N = transitionModel.getPopSize(event.type);
                    logP += Math.log(1.0/N);
                    break;

                case MIGRATE:
                    double m = transitionModel
                            .getBackwardRate(event.type, event.destType);
                    logP += Math.log(m);
                    break;

                case SAMPLE:
                    // Do nothing here: only effect of sampling event is
                    // to change the lineage counts in subsequent intervals.
                    break;
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
        Node rootNode = sbTree.getRoot(); // guarantee that rootNode is coalescent?

        // Initialise set of live nodes:
        Set<Node> liveNodes = new HashSet<>();
        liveNodes.add(rootNode);

        // Initialise lineage count:
        Integer[] lineageCount = new Integer[2];
        lineageCount[0] = ((SeedbankNodeX)rootNode).getNodeType() == 0 ? 1 : 0;
        lineageCount[1] = ((SeedbankNodeX)rootNode).getNodeType() == 1 ? 1 : 0;

        // Calculate event sequence:
        while (!liveNodes.isEmpty()) {

            SBEvent nextEvent = new SBEvent();
            nextEvent.time = Double.NEGATIVE_INFINITY;
            nextEvent.node = rootNode; // Initial assignment not significant

            // Determine next event
            for (Node node : liveNodes) {
            	if (((SeedbankNodeX)node).getNodeEvent() == NodeEvent.SAMPLE) {
            		if (node.getHeight()>nextEvent.time) {
                        nextEvent.time = node.getHeight();
                        nextEvent.kind = SBEventKind.SAMPLE;
                        nextEvent.type = ((SeedbankNodeX)node).getNodeType();
                        nextEvent.node = node;
                    }
            	} else if (((SeedbankNodeX)node).getNodeEvent() == NodeEvent.COALESCENT) {
            		if (node.getHeight()>nextEvent.time) {
                        nextEvent.time = node.getHeight();
                        nextEvent.kind = SBEventKind.COALESCE;
                        nextEvent.type = ((SeedbankNodeX)node).getNodeType();
                        nextEvent.node = node;
                    }
            	} else if (((SeedbankNodeX)node).getNodeEvent() == NodeEvent.MIGRATION) {
                    if (node.getHeight()>nextEvent.time) {
                        nextEvent.time = node.getHeight();
                        nextEvent.kind = SBEventKind.MIGRATE;
                        // destType is backward in time
                        nextEvent.destType = ((SeedbankNodeX)node).getNodeType();
                        nextEvent.type = ((SeedbankNodeX)node).getNodeToType();
                        nextEvent.node = node;
                    }
            	} else {
            		// Should not fall through.
                    throw new RuntimeException("No dummy nodes should be included here!");
            	}
            }

            // Update active node list (changeIdx) and lineage count appropriately:
            switch (nextEvent.kind) {
                case COALESCE:
                    Node leftChild = nextEvent.node.getLeft();
                    Node rightChild = nextEvent.node.getRight();
                    
                    liveNodes.remove(nextEvent.node);
                    liveNodes.add(leftChild);
                    liveNodes.add(rightChild);
                    lineageCount[nextEvent.type]++; 

                    break;

                case SAMPLE:
                	liveNodes.remove(nextEvent.node);
                    lineageCount[nextEvent.type]--;
                    break;

                case MIGRATE:
                	Node daughter = nextEvent.node.getLeft();
                	liveNodes.remove(nextEvent.node);
                    liveNodes.add(daughter);
                    
                    lineageCount[nextEvent.destType]--;
                    lineageCount[nextEvent.type]++;
                    break;
            }

            // Add event to list:
            eventList.add(nextEvent);
            lineageCountList.add(Arrays.copyOf(lineageCount, lineageCount.length));
        }

        // Reverse event and lineage count lists (order them from tips to root):
        eventList = Lists.reverse(eventList);
        lineageCountList = Lists.reverse(lineageCountList);

    }

    @Override
    public boolean requiresRecalculation() {
        return true;
    }
	
	
	// Interface requirements
	
	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
		
	}

}
