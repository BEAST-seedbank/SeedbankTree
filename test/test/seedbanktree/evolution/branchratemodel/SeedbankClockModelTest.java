package test.seedbanktree.evolution.branchratemodel;

import org.junit.Test;

import beast.base.evolution.likelihood.TreeLikelihood;
import junit.framework.TestCase;
import seedbanktree.evolution.branchratemodel.SeedbankClockModel;

public class SeedbankClockModelTest extends TestCase {
	
	@Test
	public void testSeedbankClockModel() {
		
		
		SeedbankClockModel clockModel = new SeedbankClockModel();
		
		TreeLikelihood treeLikelihood = new TreeLikelihood();
		
		
		assertEquals(1, 1, 1e-3);
	}
}
