package seedbanktree.util;

import java.io.PrintStream;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Distribution;
import seedbanktree.evolution.tree.SeedbankTree;

public class MAPTreeLogger extends Tree {
	public Input<SeedbankTree> seedbankTreeInput = new Input<>(
	        "sbTree",
	        "Seedbank tree state for logging MAP.",
	        Validate.REQUIRED);

	    public Input<Distribution> posteriorInput = new Input<>(
	        "posterior",
	        "Posterior used to identify MAP tree",
	        Validate.REQUIRED);

	    SeedbankTree currentMAPTree;
	    double maxPosterior;

	    @Override
	    public void initAndValidate() {
	        super.initAndValidate();

	        currentMAPTree = seedbankTreeInput.get().copy();
	        currentMAPTree.setTypeTrait(seedbankTreeInput.get().getTypeTrait());
	        currentMAPTree.initAndValidate();
	        maxPosterior = Double.NEGATIVE_INFINITY;
	    }

	    @Override
	    public void init(PrintStream out) {
	        currentMAPTree.init(out);
	    }

	    @Override
	    public void log(long nSample, PrintStream out) {
	        if (posteriorInput.get().getCurrentLogP()>maxPosterior) {
	            maxPosterior = posteriorInput.get().getCurrentLogP();
	            currentMAPTree.assignFrom(seedbankTreeInput.get());
	        }
	        currentMAPTree.log(nSample, out);
	    }

	    @Override
	    public void close(PrintStream out) {
	        currentMAPTree.close(out);
	    }
}
