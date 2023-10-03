package seedbanktree.evolution.tree;

import java.util.ArrayList;
import java.util.List;

import beast.base.evolution.tree.Node;

public class SeedbankNode extends Node {
	
	// Total number of changes on the branch above this node
	int nTypeChanges = 0;
    List<Integer> changeTypes = new ArrayList<Integer>();
    List<Double> changeTimes = new ArrayList<Double>();
    int nodeType = 0;
    
    // Obtain number of changes
    public int getChangeCount() {
        return nTypeChanges;
    }
    
    // Obtain destination type (reverse time) of change specified by idx
    public int getChangeType(int idx) {
        return changeTypes.get(idx);
    }
    
    // Obtain time of change specified by idx
    public double getChangeTime(int idx) {
        return changeTimes.get(idx);
    }
    
    // Get final type at end of branch above node
    public int getFinalType() {
        if (nTypeChanges>0)
            return changeTypes.get(nTypeChanges-1);
        else
            return nodeType;
    }
    
    // Get final change time
    public double getFinalChangeTime() {
        if (nTypeChanges>0)
            return changeTimes.get(nTypeChanges-1);
        else
            return getHeight(); // Node function
    }
    
    // Obtain node type
    public int getNodeType() {
        return nodeType;
    }
    
    // Set node type
    public void setNodeType(int nodeType) {
        startEditing();
        this.nodeType = nodeType;
    }
    
    // Add type change
    public void addChange(int newType, double time) {
    	// TODO: Verify that type change is valid? (actually a type change)
    	// TODO: Verify that change time added is actually in between times?
        startEditing();
        changeTypes.add(newType);
        changeTimes.add(time);
        nTypeChanges += 1;
    }
    
    // Remove type changes
    public void clearChanges() {
        startEditing();
        changeTypes.clear();
        changeTimes.clear();
        nTypeChanges = 0;
    }
    
    // Change time of type change
    public void setChangeTime(int idx, double newTime) {
        startEditing();
        changeTimes.set(idx, newTime);
    }
    
    // Change destination of type change
    // TODO: invalid? 
    public void setChangeType(int idx, int newType) {
        startEditing();
        changeTypes.set(idx, newType);
    }

    // Truncates changes
    public void truncateChanges(int newNChanges) {
        startEditing();

        while (nTypeChanges>newNChanges) {
            changeTypes.remove(nTypeChanges-1);
            changeTimes.remove(nTypeChanges-1);
            nTypeChanges -= 1;
        }
    }
    
    // Insert a new change at index idx
    public void insertChange(int idx, int newType, double newTime) {
    	// TODO: verify change?
    	// TODO: verify time of change?
        startEditing();

        if (idx>nTypeChanges)
            throw new IllegalArgumentException("Index to insertChange() out of range.");

        changeTimes.add(idx, newTime);
        changeTypes.add(idx, newType);
        nTypeChanges += 1;
    }
    
    // Remove change
    public void removeChange(int idx) {
        startEditing();

        if (idx>=nTypeChanges)
            throw new IllegalArgumentException("Index to removeChange() out of range.");

        changeTimes.remove(idx);
        changeTypes.remove(idx);
        nTypeChanges -= 1;

    }
    
    /**
     * @return shallow copy of node
     */
    public SeedbankNode shallowCopy() {
    	SeedbankNode node = new SeedbankNode();
        node.height = height;
        node.parent = parent;        
        node.children.addAll(children);

        node.nTypeChanges = nTypeChanges;
        node.changeTimes.addAll(changeTimes);
        node.changeTypes.addAll(changeTypes);
        node.nodeType = nodeType;
                
        node.labelNr = labelNr;
        node.metaDataString = metaDataString;
        node.ID = ID;

        return node;
    }
    
    /**
     * **************************
     * Methods ported from Node *
     ***************************
     */


    /**
     * @return (deep) copy of node
     */
    @Override
    public SeedbankNode copy() {
    	SeedbankNode node = new SeedbankNode();
        node.height = height;
        node.labelNr = labelNr;
        node.metaDataString = metaDataString;
        node.parent = null;
        node.ID = ID;
        node.nTypeChanges = nTypeChanges;
        node.changeTimes.addAll(changeTimes);
        node.changeTypes.addAll(changeTypes);
        node.nodeType = nodeType;
        if (getLeft()!=null) {
            node.setLeft(getLeft().copy());
            node.getLeft().setParent(node);
            if (getRight()!=null) {
                node.setRight(getRight().copy());
                node.getRight().setParent(node);
            }
        }
        return node;
    }

    /**
     * assign values from a tree in array representation *
     * @param nodes
     * @param node
     */
    @Override
    public void assignFrom(Node[] nodes, Node node) {
        height = node.getHeight();
        labelNr = node.getNr();
        metaDataString = node.metaDataString;
        parent = null;
        ID = node.getID();
        
        SeedbankNode mtNode = (SeedbankNode)node;
        nTypeChanges = mtNode.nTypeChanges;
        changeTimes.clear();
        changeTimes.addAll(mtNode.changeTimes);
        changeTypes.clear();
        changeTypes.addAll(mtNode.changeTypes);
        nodeType = mtNode.nodeType;
        
        if (node.getLeft()!=null) {
            setLeft(nodes[node.getLeft().getNr()]);
            getLeft().assignFrom(nodes, node.getLeft());
            getLeft().setParent(this);
            if (node.getRight()!=null) {
                setRight(nodes[node.getRight().getNr()]);
                getRight().assignFrom(nodes, node.getRight());
                getRight().setParent(this);
            }
        }
    }
}
