package seedbanktree.evolution.tree;

import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;

@Description("Class to initialize a SeedbankTree from single child newick tree with type metadata")
public class SeedbankTreeFromNewick extends SeedbankTree implements StateNodeInitialiser {

	public Input<String> newickStringInput = new Input<>("value",
            "Tree in Newick format.", Validate.REQUIRED);

    public Input<Boolean> adjustTipHeightsInput = new Input<>("adjustTipHeights",
            "Adjust tip heights in tree? (Default false).", false);
    
    public Input <Alignment> dataInput = new Input<> ("data", 
    		"Specifies the sequences represented by the leaves in the tree.");
    
    
	public SeedbankTreeFromNewick() {
    	m_initial.setRule(Validate.REQUIRED);
    }
	
	@Override
	public void initAndValidate() {
		// seedbankTree inputs
    	typeLabel = typeLabelInput.get();
        activeTypeName = activeTypeNameInput.get();
        dormantTypeName = dormantTypeNameInput.get();
        
//        processTraits(m_traitList.get());
        
        // parse newick string
        TreeParser parser = new TreeParser();
        parser.initByName(
                "IsLabelledNewick", true,
                "adjustTipHeights", adjustTipHeightsInput.get(),
                "singlechild", true,
                "newick", newickStringInput.get(),
                "taxa", dataInput.get());

        initFromFlatTree(parser, true);
	}
	
    // Methods for StateNodeInitialiser interface

	@Override
    public void initStateNodes() {
        assert(m_initial != null);
        if (!(m_initial.get() instanceof SeedbankTree)) {
        	throw new IllegalArgumentException("Attempted to use "
                    + "seedbank tree initialiser on regular tree object");
        }
    	m_initial.get().assignFromWithoutID(this);
    }

    @Override
    public void getInitialisedStateNodes(final List<StateNode> stateNodes) {
    	assert(m_initial != null);
        if (!(m_initial.get() instanceof SeedbankTree)) {
        	throw new IllegalArgumentException("Attempted to use "
                    + "seedbank tree initialiser on regular tree object");
        }
        stateNodes.add(m_initial.get());
    }

}
