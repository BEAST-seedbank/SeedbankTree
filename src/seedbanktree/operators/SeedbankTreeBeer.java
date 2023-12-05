package seedbanktree.operators;

import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.inference.Operator;
import beast.base.util.Randomizer;
import multitypetree.operators.BeerliFelsenstein.Event;
import seedbanktree.evolution.tree.SeedbankTree;
import seedbanktree.evolution.tree.SeedbankNode;

public class SeedbankTreeBeer extends SeedbankTreeOperator {
	
	private enum EventType {MIGRATION, COALESCENCE, SAMPLE};
    private class Event {
        EventType eventType;
        double time;
        private EventType type;
        Node node;
        int thisDeme, prevDeme; //prevDeme only for migration events
    }

    @Override
    public void initAndValidate() {
        sbTree = seedbankTreeInput.get();
//            trModel = transitionModelInput.get();
    }

	@Override
	public double proposal() {
		double logHR = 0.0;
		
		// Select non-root node at random
        Node node = sbTree.getNode(Randomizer.nextInt(sbTree.getNodeCount())-1);
        SeedbankNode sbNode = (SeedbankNode)node;
        
        // Keep copies of node and its sister for reverse move prob calculation
        SeedbankNode sbNodeOld = sbNode.shallowCopy();
        SeedbankNode sbNodeOldSis = (SeedbankNode)getOtherChild(node.getParent(), node);
        
        // Assemble partial event list
        List<Event> eventList = getPartialEventList(node);
        double oldRootHeight = sbTree.getRoot().getHeight();
        
        // Topology changes to turn tree into partial tree
        Node nodeParent = node.getParent();
        if (!nodeParent.isRoot())
            disconnectBranch(node);
        else { //if node's parent is root
            Node sister = getOtherChild(nodeParent, node);
            nodeParent.removeChild(sister);
        }
        
        // Pre-calculate total lineage migration propensities
        double [] props = new double[2];
        for (int d=0; d<2; d++) {
        	props[d] = 0.0;
            for (int dp=0; dp<2; dp++) {
                if (d==dp)
                    continue;
                props[d] += trModel.getBackwardRate(d, dp);
            }
        }
		
        List<Set<Node>> nodesOfType = Lists.newArrayList();
        nodesOfType.add(new HashSet<Node>());
        nodesOfType.add(new HashSet<Node>());

        sbNode.clearChanges();
        
        // simulation starts
        double coalTime = Double.NaN;
        // iterate through partial event list
        // (simulation beneath the root)
        for (int eidx=0; eidx<eventList.size(); eidx++) { 
            Event event = eventList.get(eidx);
            double intervalEndTime;
            if (eidx<eventList.size()-1) // if not last event before root
                intervalEndTime = eventList.get(eidx+1).time;
            else
                intervalEndTime = oldRootHeight;
            
            switch (event.type) {
                case COALESCENCE:
//                    nodesOfType.get(event.thisDeme).removeAll(event.node.getChildren());
//                    nodesOfType.get(event.thisDeme).add(event.node);
                	nodesOfType.get(1).removeAll(event.node.getChildren()); //coalescent -> must be active = 1
                    nodesOfType.get(1).add(event.node);
                    break;
                    
                case SAMPLE:
                    nodesOfType.get(event.thisDeme).add(event.node);
                    break;
                    
                case MIGRATION:
                    nodesOfType.get(event.prevDeme).remove(event.node);
                    nodesOfType.get(event.thisDeme).add(event.node);
                    break;
            }
            
            double t = Math.max(event.time, node.getHeight());
            
            // Early exit 
            if (t >= intervalEndTime)
                continue; // skip until t < intervalEndTime
            	// i.e. Math.max(event.time, node.getHeight()) < intervalEndTime
            	// i.e. intervalEndTime > event.time AND intervalEndTime > node.getHeight()
            	// Iterate until point in eventList where [event, event+1] contains node

            int deme = sbNode.getNodeType();
            
            while (true) { // simulate until coalesce
                
                // Calculate coalescent propensity
            	// only if Active
            	double coalProp;
            	if (deme == 1)
            		coalProp = nodesOfType.get(1).size()/trModel.getPopSize(1); 
            	else 
            		coalProp = 0;
            	
            	
                // Select event time
                double dt = Randomizer.nextExponential(coalProp + props[deme]);
                t += dt;
                if (t > intervalEndTime) // if time of next event drawn is greater than time of event+1
                    break;
            
                // HR waiting time contribution
                // TODO: waiting time propensity change?
                logHR += -(coalProp + props[deme])*dt;
                
                double u = Randomizer.nextDouble()*(coalProp + props[deme]);
                if (u<coalProp) {
                    // Coalescence

                    // Select edge to coalesce with
                    Node coalNode = (Node)selectRandomElement(nodesOfType.get(deme));
                    
                    // HR event contribution
                    logHR += Math.log(1.0/trModel.getPopSize(deme));
                    
                    // Implement coalescence
                    coalTime = t;                    
                    connectBranch(node, coalNode, coalTime);

                    break;
                
                } else {
                    // Migration
                
                    u -= coalProp;
                    int toDeme;
                    for (toDeme = 0; toDeme<2; toDeme++) {
                        if (toDeme == deme)
                            continue;
                    
                        u -= trModel.getBackwardRate(deme, toDeme);
                        if (u<0)
                            break;
                    }
                
                    // HR event contribution
                    logHR += Math.log(trModel.getBackwardRate(deme, toDeme));

                    // Implelent migration
                    sbNode.addChange(toDeme, t);
                    deme = toDeme;
                }
            }
            
            // Continue to next interval if no coalescence has occurred
            if (!Double.isNaN(coalTime))
                break;
        }
        
        if (Double.isNaN(coalTime)) {
            
            // Continue simulation beyond old root of tree
            double t = oldRootHeight;
            
            int deme = sbNode.getFinalType();
            SeedbankNode sbNodeSis = (SeedbankNode)eventList.get(eventList.size()-1).node; //oldroot?
            int demeSis = sbNodeSis.getFinalType();
            
            while (true) {
                
                // Calculate coalescent propensity
            	// TODO: change 
                double coalProp;
                if (deme == demeSis && deme == 1)
                    coalProp = 1.0/trModel.getPopSize(deme);
                else
                    coalProp = 0.0;
                
                double totalProp = coalProp + props[deme] + props[demeSis]; //TODO: not sure
                double dt = Randomizer.nextExponential(totalProp);
                
                // HR waiting time contribution
                logHR += -totalProp*dt;
                
                t += dt;
                
                double u = Randomizer.nextDouble()*totalProp;
                
                if (u <coalProp) { 
                    // Coalescence
                    
                    logHR += Math.log(1.0/trModel.getPopSize(deme));
                    
                    coalTime = t;
                    nodeParent.addChild(sbNodeSis);
                    sbTree.setRoot(nodeParent);
                    
                    break; // end simulation, coalescent with last node completed
                    
                } else {
                    // Migration
                    
                    u -= coalProp;
                    
                    if (u<props[deme]) {
                        // Migration in main lineage
                        
                        int toDeme;
                        for (toDeme=0; toDeme<2; toDeme++) {
                            if (toDeme == deme)
                                continue;
                            
                            u -= trModel.getBackwardRate(deme, toDeme);
                            if (u<0)
                                break;
                        }
                        
                        // HR contribution
                        logHR += Math.log(trModel.getBackwardRate(deme, toDeme));
                        
                        sbNode.addChange(toDeme, t);
                        deme = toDeme;
                    } else {
                        // Migration in sister lineage
                        
                        int toDeme;
                        for (toDeme=0; toDeme<2; toDeme++) {
                            if (toDeme == demeSis)
                                continue;
                            
                            u -= trModel.getBackwardRate(demeSis, toDeme);
                            if (u<0)
                                break;
                        }
                        
                        // HR contribution
                        logHR += Math.log(trModel.getBackwardRate(demeSis, toDeme));
                        
                        sbNodeSis.addChange(toDeme, t);
                        demeSis = toDeme;
                    }
                }
            }
        }
        
		return logHR;
	}
	
    /**
     * Assemble and return list of events excluding those on the edge between
     * node and its parent.
     * 
     * @param excludedNode Tree node indicating edge to exclude.
     * @return event list
     */
    private List<Event> getPartialEventList(Node excludedNode) {
        List<Event> eventList = Lists.newArrayList();
        
        
        // Collect all events
        for (Node node : sbTree.getNodesAsArray()) {
            
            if (node == excludedNode)
                continue;

            SeedbankNode sbNode = (SeedbankNode)node;
            
            Event event = new Event();
            event.time = node.getHeight();
            event.node = node;
            event.thisDeme = sbNode.getNodeType();
            if (node.isLeaf())
                event.type = EventType.SAMPLE;
            else {
                if (!node.getChildren().contains(excludedNode))
                    event.type = EventType.COALESCENCE;
            }
            eventList.add(event);
            
            int thisDeme = sbNode.getNodeType();
            int prevDeme; //prevDeme only for migration events
            for (int i=0; i<sbNode.getChangeCount(); i++) {
                prevDeme = thisDeme;
                thisDeme = sbNode.getChangeType(i);
                
                Event changeEvent = new Event();
                changeEvent.type = EventType.MIGRATION;
                changeEvent.time = sbNode.getChangeTime(i);
                changeEvent.node = node;
                changeEvent.thisDeme = thisDeme;
                changeEvent.prevDeme = prevDeme; 
                eventList.add(changeEvent);
            }
            
        }
        
        // Sort events according to times
        // from smallest time to biggest time
        // leaves to root
        Collections.sort(eventList, new Comparator<Event>() {

            @Override
            public int compare(Event e1, Event e2) {
                if (e1.time < e2.time)
                    return -1;
                
                if (e1.time > e2.time)
                    return 1;
                
                return 0;
            }
        });

        return eventList;
    }
    
    /**
     * Select element at random from set.
     * 
     * @param set
     * @return Object
     */
    public Object selectRandomElement(Set set) {
        return set.toArray()[Randomizer.nextInt(set.size())];
    }

}
