package seedbanktree.evolution.tree;

import java.util.ArrayList;
import java.util.List;

import beast.base.evolution.tree.Node;

public class SeedbankNodeX extends Node {
	
	int nodeType;
//    int nodeFromType;
	
	// The type post migration, time forward. 
    int nodeToType;
    
    // 0: Sample, 1: Migration, 2: Coalescent, 3: Empty (Migration dummy node)
    public enum NodeEvent {
    	SAMPLE, MIGRATION, COALESCENT, DUMMY;
    }
    
    NodeEvent nodeEvent;
    
    // Obtain node type
    public int getNodeType() {
        return nodeType;
    }
    
//    public int getNodeFromType() {
//        return nodeFromType;
//    }
    
    public int getNodeToType() {
        return nodeToType;
    }
    
    public NodeEvent getNodeEvent() {
    	return nodeEvent;
    }
    
    public void setNodeType(int nodeType) {
        startEditing();
        this.nodeType = nodeType;
    }
    
//    public void setNodeFromType(int nodeFromType) {
//        startEditing();
//        this.nodeFromType = nodeFromType;
//    }
    
    public void setNodeToType(int nodeToType) {
        startEditing();
        this.nodeToType = nodeToType;
    }
    
    public void setNodeEvent(NodeEvent nodeEvent) {
    	startEditing();
    	this.nodeEvent = nodeEvent;
    }
 

    /**
     * @return shallow copy of node
     */
    public SeedbankNodeX shallowCopy() {
    	SeedbankNodeX node = new SeedbankNodeX();
        node.height = height;
        node.parent = parent;        
        node.children.addAll(children);

        node.nodeType = nodeType;
//        node.nodeFromType = nodeFromType;
        node.nodeToType = nodeToType;
        node.nodeEvent = nodeEvent;
                
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
    public SeedbankNodeX copy() {
    	SeedbankNodeX node = new SeedbankNodeX();
        node.height = height;
        node.labelNr = labelNr;
        node.metaDataString = metaDataString;
        node.parent = null;
        node.ID = ID;
        node.nodeType = nodeType;
//        node.nodeFromType = nodeFromType;
        node.nodeToType = nodeToType;
        node.nodeEvent = nodeEvent;
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
        
        SeedbankNodeX mtNode = (SeedbankNodeX)node;
        nodeType = mtNode.nodeType;
//        nodeFromType = mtNode.nodeFromType;
        nodeToType = mtNode.nodeToType;
        nodeEvent = mtNode.nodeEvent;
        
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
