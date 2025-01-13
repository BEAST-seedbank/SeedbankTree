package seedbanktree.util;

import java.io.PrintStream;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import seedbanktree.evolution.tree.SeedbankNode;
import seedbanktree.evolution.tree.SeedbankTree;

public class DormancyTreeLogger extends BEASTObject implements Loggable {

    public Input<SeedbankTree> seedbankTreeInput = new Input<>(
            "sbTree",
            "Seedbanktree tree to log.",
            Input.Validate.REQUIRED);

    SeedbankTree sbTree;

    @Override
    public void initAndValidate() {
        sbTree = seedbankTreeInput.get();
    }

    @Override
    public void init(PrintStream out) {
        sbTree.init(out);
    }

    @Override
    public void log(long nSample, PrintStream out) {
    	
        // Set up metadata string
        for (Node node : sbTree.getNodesAsArray()) {
            SeedbankNode sbNode = (SeedbankNode)node;
            
            double dormantPercentage = 0.0;
            
            if (!sbNode.isRoot()) {
            	double dormantLength = 0.0;

                int thisType = sbNode.getNodeType();
                double lastTime = sbNode.getHeight();
                for (int i = 0; i < sbNode.getChangeCount(); i++) {
                    double nextTime = sbNode.getChangeTime(i);
                    dormantLength += thisType == 0 ? (nextTime - lastTime) : 0;
                    lastTime = nextTime;
                    thisType = 1 - thisType;
                }
                   
                dormantPercentage = dormantLength / sbNode.getLength();
            }
            
            sbNode.metaDataString = String.format("lambda=%.3f", dormantPercentage);
        }
            
        out.print("tree STATE_" + nSample + " = ");
        out.print(sbTree.getRoot().toSortedNewick(new int[1], true));
        out.print(";");
    }

    @Override
    public void close(PrintStream out) {
        sbTree.close(out);
    }
}
