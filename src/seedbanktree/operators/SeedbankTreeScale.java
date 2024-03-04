package seedbanktree.operators;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.inference.Operator;
import beast.base.util.Randomizer;
import seedbanktree.evolution.tree.SeedbankNode;
import seedbanktree.evolution.tree.SeedbankTree;
import seedbanktree.evolution.tree.TransitionModel;

public class SeedbankTreeScale extends Operator{

    public Input<SeedbankTree> seedbankTreeInput = new Input<>(
            "seedbankTree", "Seedbank tree on which to operate.",
            Validate.REQUIRED);
            
//    public Input<TransitionModel> transitionModelInput = new Input<>(
//        "transitionModel",
//        "Transition model for proposal distribution",
//        Input.Validate.REQUIRED);
        
    public Input<Double> scaleFactorInput = new Input<>("scaleFactor",
    		"Scaling is restricted to the range [1/scaleFactor, scaleFactor]");
        
    protected SeedbankTree sbTree;
//    protected TransitionModel trModel;

    @Override
    public void initAndValidate() {
    	sbTree = seedbankTreeInput.get();
//        trModel = transitionModelInput.get();
    }

	@Override
	public double proposal() {
		// Choose scale factor:
        double u = Randomizer.nextDouble();
        double f = u*scaleFactorInput.get()+(1.0-u)/scaleFactorInput.get();

        // Keep track of Hastings ratio:
        double logf = Math.log(f);
        double logHR = -2*logf;
       
        // Scale colour change times on external branches:        
        for (Node leaf : sbTree.getExternalNodes()) {
            for (int c = 0; c<((SeedbankNode)leaf).getChangeCount(); c++) {
                double oldTime = ((SeedbankNode)leaf).getChangeTime(c);
                ((SeedbankNode)leaf).setChangeTime(c, f*oldTime);
                logHR += logf;
            }
        }

        // Scale internal node heights and colour change times:
        for (Node node : sbTree.getInternalNodes()) {
            
            node.setHeight(node.getHeight()*f);
            logHR += logf;
            
            for (int c = 0; c<((SeedbankNode)node).getChangeCount(); c++) {
                double oldTime = ((SeedbankNode)node).getChangeTime(c);
                ((SeedbankNode)node).setChangeTime(c, f*oldTime);
                logHR += logf;
            }
        }
        
        // Reject invalid tree scalings:
        if (f<1.0) {
            for (Node leaf : sbTree.getExternalNodes()) {
                if (leaf.getParent().getHeight()<leaf.getHeight())
                    return Double.NEGATIVE_INFINITY;
                
                if (((SeedbankNode)leaf).getChangeCount()>0
                        && ((SeedbankNode)leaf).getChangeTime(0)<leaf.getHeight())
                        return Double.NEGATIVE_INFINITY;
            }
        }
         
        return logHR;
	}

}
