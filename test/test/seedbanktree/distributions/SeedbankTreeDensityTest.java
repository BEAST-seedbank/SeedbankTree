package test.seedbanktree.distributions;

import org.junit.Test;

import beast.base.evolution.tree.Node;
import beast.base.inference.parameter.RealParameter;
import junit.framework.TestCase;
import seedbanktree.distributions.SeedbankTreeDensity;
import seedbanktree.evolution.tree.SeedbankTree;
import seedbanktree.evolution.tree.SeedbankTreeFromNewick;
import seedbanktree.evolution.tree.TransitionModel;

public class SeedbankTreeDensityTest extends TestCase {

	/**
	 * Test #1 of calculateLogP method of class SeedbankTreeDensity.
	 */
	@Test
	public void testCalculateLogPOne() {
		// Assemble test tree:
		String newickStr =
                        "(((((A[&state=1]:12.0,B[&state=1]:12.0)[&state=1]:41.0,((C[&state=1]:10.0)"
                        + "[&state=0]:17.0)[&state=1]:26.0)[&state=1]:10.0)[&state=0]:15.0)[&state=1]:3.0,"
                        + "((D[&state=1]:23.0)[&state=0]:10.0)[&state=1]:48.0)[&state=1]:0;";

		SeedbankTreeFromNewick sbTreeInitializer = new SeedbankTreeFromNewick();
		SeedbankTree sbTree = new SeedbankTree();
		
		sbTreeInitializer.initByName(
                        "value", newickStr,
                        "typeLabel", "state",
                        "initial", sbTree);
		
		sbTreeInitializer.initStateNodes();

		// Assemble transition model:
		RealParameter rate = new RealParameter();
		rate.initByName("value","1.0");
		RealParameter K = new RealParameter();
		K.initByName("value","1.0");
		RealParameter activeSize = new RealParameter();
		activeSize.initByName(
						"value", "100.0",
						"estimate", false);
		
		TransitionModel transitionModel = new TransitionModel();
		transitionModel.initByName(
                        "rate", rate,
                        "K", K,
                        "activeSize", activeSize);

		// Set up likelihood instance:
		SeedbankTreeDensity likelihood = new SeedbankTreeDensity();
		likelihood.initByName(
                        "transitionModel", transitionModel,
                        "seedbankTree", sbTree,
                        "checkValidity", true);

		double expResult = -243.675; // Calculated by hand
		double result = likelihood.calculateLogP();
		
		System.out.println(result);
		assertEquals(expResult, result, 1e-3);
		
	}
	
	/**
	 * Test #2 of calculateLogP method of class SeedbankTreeDensity.
	 */
	@Test
	public void testCalculateLogPTwo() {
		String newickStr =
                        "((((A[&state=1]:31.0,B[&state=1]:31.0)[&state=1]:5.0)[&state=0]:15.0)"
                        + "[&state=1]:3.0,(((C[&state=1]:8.0)[&state=0]:7.0)[&state=1]:11.0,"
                        + "D[&state=1]:7.0)[&state=1]:33.0)[&state=1]:0.0;";

		SeedbankTreeFromNewick sbTreeInitializer = new SeedbankTreeFromNewick();
		SeedbankTree sbTree = new SeedbankTree();
		
		sbTreeInitializer.initByName(
                        "value", newickStr,
                        "typeLabel", "state",
                        "initial", sbTree);
		
		sbTreeInitializer.initStateNodes();

		// Assemble transition model:
		RealParameter rate = new RealParameter();
		rate.initByName("value","2.0");
		RealParameter K = new RealParameter();
		K.initByName("value","1.25");
		RealParameter activeSize = new RealParameter();
		activeSize.initByName(
						"value", "90.0",
						"estimate", false);
		
		TransitionModel transitionModel = new TransitionModel();
		transitionModel.initByName(
                        "rate", rate,
                        "K", K,
                        "activeSize", activeSize);

		// Set up likelihood instance:
		SeedbankTreeDensity likelihood = new SeedbankTreeDensity();
		likelihood.initByName(
                        "transitionModel", transitionModel,
                        "seedbankTree", sbTree,
                        "checkValidity", true);

		double expResult = -325.96; // Calculated by hand
		double result = likelihood.calculateLogP();
		
		System.out.println(result);
		assertEquals(expResult, result, 1e-3);
		
//		for (Node n : sbTree.getNodesAsArray()) {
//			System.out.println(n.getID());
//			System.out.println(n.getHeight());
//		}
		
	}
}
