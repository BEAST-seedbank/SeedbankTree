package test.seedbanktree.evolution.tree;

import org.junit.Ignore;
import org.junit.Test;

import beast.base.evolution.tree.Node;
import junit.framework.TestCase;
import seedbanktree.evolution.tree.SeedbankNode;
import seedbanktree.evolution.tree.SeedbankTree;
import seedbanktree.evolution.tree.SeedbankTreeFromNewick;

public class IsValidTest extends TestCase {
	
	@Test
	public void testCheckValidity() {
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
		
		assertTrue("This tree is valid", sbTree.isValid());
	}
	
	@Test
	public void testCheckValidityTwo() {
		// Assemble test tree:
		String newickStr =
                "((A[&state=0]:10)[&state=1]:10, B[&state=1]:10)[&state=1]:0;";

		SeedbankTreeFromNewick sbTreeInitializer = new SeedbankTreeFromNewick();
		SeedbankTree sbTree = new SeedbankTree();
		
		sbTreeInitializer.initByName(
                        "value", newickStr,
                        "typeLabel", "state",
                        "initial", sbTree);
		
		sbTreeInitializer.initStateNodes();
		
		assertTrue("This tree is valid", sbTree.isValid());
	}
	
	@Test
	public void testCheckValidityThree() {
		// Assemble test tree:
		String newickStr =
                "((A[&state=0]:10)[&state=0]:10, B[&state=1]:10)[&state=1]:0;";

		SeedbankTreeFromNewick sbTreeInitializer = new SeedbankTreeFromNewick();
		SeedbankTree sbTree = new SeedbankTree();
		
		sbTreeInitializer.initByName(
                        "value", newickStr,
                        "typeLabel", "state",
                        "initial", sbTree);
		
		sbTreeInitializer.initStateNodes();
		
		assertFalse("This tree has an invalid migration", sbTree.isValid());
	}
	
	@Test
	public void testCheckValidityFour() {
		// Assemble test tree:
		String newickStr =
                "((A[&state=1]:10)[&state=0]:10, B[&state=1]:10)[&state=1]:0;";

		SeedbankTreeFromNewick sbTreeInitializer = new SeedbankTreeFromNewick();
		SeedbankTree sbTree = new SeedbankTree();
		
		sbTreeInitializer.initByName(
                        "value", newickStr,
                        "typeLabel", "state",
                        "initial", sbTree);
		
		sbTreeInitializer.initStateNodes();
		
		assertFalse("This tree has an invalid coalescent", sbTree.isValid());
	}
	
	@Test
	public void testCheckValidityFive() {
		// Assemble test tree:
		String newickStr =
                "((A[&state=1]:10)[&state=1]:10, B[&state=1]:10)[&state=1]:0;";

		SeedbankTreeFromNewick sbTreeInitializer = new SeedbankTreeFromNewick();
		SeedbankTree sbTree = new SeedbankTree();
		
		sbTreeInitializer.initByName(
                        "value", newickStr,
                        "typeLabel", "state",
                        "initial", sbTree);
		
		sbTreeInitializer.initStateNodes();
		
		assertFalse("This tree has an invalid migration", sbTree.isValid());
	}
	
	@Test
	public void testCheckValiditySix() {
		// Assemble test tree:
		String newickStr =
                "((A[&state=1]:10)[&state=0]:10, B[&state=0]:10)[&state=0]:0;";

		SeedbankTreeFromNewick sbTreeInitializer = new SeedbankTreeFromNewick();
		SeedbankTree sbTree = new SeedbankTree();
		
		sbTreeInitializer.initByName(
                        "value", newickStr,
                        "typeLabel", "state",
                        "initial", sbTree);
		
		sbTreeInitializer.initStateNodes();
		
		assertFalse("This tree has an invalid coalescent", sbTree.isValid());
	}
	
	@Test
	public void testCheckValiditySeven() {
		// Assemble test tree:
		String newickStr =
                        "(((A[&state=0]:10)[&state=0]:10)[&state=1]:10, B[&state=1]:10)[&state=1]:0;";

		SeedbankTreeFromNewick sbTreeInitializer = new SeedbankTreeFromNewick();
		SeedbankTree sbTree = new SeedbankTree();
		
		sbTreeInitializer.initByName(
                        "value", newickStr,
                        "typeLabel", "state",
                        "initial", sbTree);
		
		sbTreeInitializer.initStateNodes();
		
		assertFalse("This tree has an invalid migration", sbTree.isValid());
	}
	
	@Test
	public void testCheckValidityEight() {
		// Assemble test tree:
		String newickStr =
                        "(((A[&state=1]:10)[&state=0]:10)[&state=1]:-1, B[&state=1]:10)[&state=1]:0;";

		SeedbankTreeFromNewick sbTreeInitializer = new SeedbankTreeFromNewick();
		SeedbankTree sbTree = new SeedbankTree();
		
		sbTreeInitializer.initByName(
                        "value", newickStr,
                        "typeLabel", "state",
                        "initial", sbTree);
		
		sbTreeInitializer.initStateNodes();
		
		assertFalse("This tree has invalid heights", sbTree.isValid());
	}
	
	@Test
	public void testCheckValidityNine() {
		// Assemble test tree:
		String newickStr =
                        "(((A[&state=1]:10)[&state=0]:-1)[&state=1]:10, B[&state=1]:10)[&state=1]:0;";

		SeedbankTreeFromNewick sbTreeInitializer = new SeedbankTreeFromNewick();
		SeedbankTree sbTree = new SeedbankTree();
		
		sbTreeInitializer.initByName(
                        "value", newickStr,
                        "typeLabel", "state",
                        "initial", sbTree);
		
		sbTreeInitializer.initStateNodes();
		
		assertFalse("This tree has invalid heights", sbTree.isValid());
	}
	
	@Test
	public void testCheckValidityTen() {
		// Assemble test tree:
		String newickStr =
                        "((((A[&state=1]:10)[&state=0]:10)[&state=1]:10, B[&state=1]:10)[&state=1]:-1,C[&state=1]:10)[&state=1]:0;";

		SeedbankTreeFromNewick sbTreeInitializer = new SeedbankTreeFromNewick();
		SeedbankTree sbTree = new SeedbankTree();
		
		sbTreeInitializer.initByName(
                        "value", newickStr,
                        "typeLabel", "state",
                        "initial", sbTree);
		
		sbTreeInitializer.initStateNodes();
		
//		for (Node n : sbTree.getNodesAsArray()) {
//			SeedbankNode sbNode = (SeedbankNode) n;
//			System.out.println(String.format("ID: %s, type: %d, height: %f", sbNode.getID(), sbNode.getNodeType(), sbNode.getHeight()));
//			System.out.println(sbNode.changeTypes);
//			System.out.println(sbNode.changeTimes);
//		}
//		System.out.println(sbTree.isValid());
		
		assertFalse("This tree has invalid heights", sbTree.isValid());
	}

}


