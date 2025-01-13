package seedbanktree.util;

import java.io.PrintStream;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import seedbanktree.evolution.tree.SeedbankNode;
import seedbanktree.evolution.tree.SeedbankTree;

public class TypedNodeTreeLogger extends BEASTObject implements Loggable {

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
            String type = (sbNode.getNodeType() == 0) ? "dormant" : "active";
            sbNode.metaDataString = sbTree.getTypeLabel()
                    + "=\""
                    + type
                    + "\"";
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
