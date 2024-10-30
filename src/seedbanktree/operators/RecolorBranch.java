package seedbanktree.operators;

import beast.base.core.Description;
import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import seedbanktree.evolution.tree.SeedbankNode;
import seedbanktree.evolution.tree.SeedbankTree;

@Description("Proposes a new coloring for a randomly selected branch.")
public class RecolorBranch extends UniformizationRetypeOperator {

    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
	@Override
    public double proposal() {
		Node node;
		do {
            node = sbTree.getNode(Randomizer.nextInt(sbTree.getNodeCount()));
        } while (node.isRoot());
		
		
		double logHR = 0.0;
		logHR += getBranchTypeProb(node);
		try {
			logHR -= retypeBranch(node);
		} catch (NoValidPathException e) {
			return Double.NEGATIVE_INFINITY;
		}
		
		return logHR;
    }

}