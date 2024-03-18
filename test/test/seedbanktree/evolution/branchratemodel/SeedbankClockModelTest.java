package test.seedbanktree.evolution.branchratemodel;

import org.junit.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.JukesCantor;
import beast.base.inference.parameter.RealParameter;
import junit.framework.TestCase;
import seedbanktree.evolution.branchratemodel.SeedbankClockModel;
import seedbanktree.evolution.tree.SeedbankTree;
import seedbanktree.evolution.tree.SeedbankTreeFromNewick;

public class SeedbankClockModelTest extends TestCase {
	
	@Test
	public void testSeedbankClockModel() {
		// Assemble data:
		Sequence sample0 = new Sequence("sample0", "AATCGGAGTT");
		Sequence sample1 = new Sequence("sample1", "ATTACCTCAT");
		Sequence sample2 = new Sequence("sample2", "CACTGTCATC");
		
		Alignment data = new Alignment();
		data.initByName(
				"sequence", sample0,
				"sequence", sample1,
				"sequence", sample2,
				"dataType", "nucleotide");
		
		// Assemble test tree:
		String newickStr =
                        "((sample0[&state=1]:1.0,((sample1[&state=1]:1.0)[&state=0]:1.0)"
                        + "[&state=1]:1.0)[&state=1]:1.0,((sample2[&state=1]:1.0)"
                        + "[&state=0]:1.0)[&state=1]:1.0)[&state=1]:0.0;";
		
		SeedbankTreeFromNewick sbTreeInitializer = new SeedbankTreeFromNewick();
		SeedbankTree sbTree = new SeedbankTree();
		
		sbTreeInitializer.initByName(
                        "value", newickStr,
                        "typeLabel", "state",
                        "initial", sbTree,
                        "data", data);
		
		sbTreeInitializer.initStateNodes();
		
		// Assemble clock model
		RealParameter activeRate = new RealParameter();
		activeRate.initByName("value", "1.0", "estimate", false);
		
		RealParameter dormantRate = new RealParameter();
		dormantRate.initByName("value", "0.5", "estimate", false);
		
		SeedbankClockModel clockModel = new SeedbankClockModel();
		clockModel.initByName(
						"activeRate", activeRate,
						"dormantRate", dormantRate,
						"tree", sbTree);
		
		JukesCantor substModel = new JukesCantor();
		substModel.initAndValidate();
	
		SiteModel siteModel = new SiteModel();
		siteModel.initByName("substModel", substModel);
		
		TreeLikelihood treeLikelihood = new TreeLikelihood();
		treeLikelihood.initByName(
						"data", data,
						"tree", sbTree,
						"siteModel", siteModel,
						"branchRateModel", clockModel,
						"implementation", "beast.evolution.likelihood.TreeLikelihood");
		
		
		double expResult = 1; // Calculated by hand
		double result = treeLikelihood.calculateLogP();
		
		System.out.println(result);
		assertEquals(expResult, result, 1e-3);
	}
}
