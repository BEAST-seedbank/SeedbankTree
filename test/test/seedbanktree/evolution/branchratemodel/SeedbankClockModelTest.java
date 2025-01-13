package test.seedbanktree.evolution.branchratemodel;

import static org.junit.Assert.assertEquals;

import org.junit.Ignore;
import org.junit.Test;

import beast.base.evolution.tree.Node;
import beast.base.core.Function.Constant;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.JukesCantor;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.GTR;
import beast.base.inference.parameter.RealParameter;
import seedbanktree.evolution.branchratemodel.SeedbankClockModel;
import seedbanktree.evolution.tree.SeedbankNode;
import seedbanktree.evolution.tree.SeedbankTree;
import seedbanktree.evolution.tree.SeedbankTreeFromNewick;

public class SeedbankClockModelTest {
	
	@Test
	public void testSeedbankClockModel() {
		// Single site
		// Assemble data:
		Sequence sample0 = new Sequence("sample0", "A");
		Sequence sample1 = new Sequence("sample1", "C");
		Sequence sample2 = new Sequence("sample2", "G");
		
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
		RealParameter dormantScaling = new RealParameter();
		dormantScaling.initByName("value", "0.5");
		
		SeedbankClockModel clockModel = new SeedbankClockModel();
		clockModel.initByName(
						"dormantScaling", dormantScaling,
						"tree", sbTree);
		
		// Assemble substitution model, site model
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
				
		double expResult = -4.17099708587562; // Calculated by hand
		double result = treeLikelihood.calculateLogP();
		
		System.out.println(result);
		assertEquals(expResult, result, 1e-5);
	}
	
	@Test
	public void testSeedbankClockModel2() {
		// 10 sites
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
		RealParameter dormantScaling = new RealParameter();
		dormantScaling.initByName("value", "0.5");
		
		SeedbankClockModel clockModel = new SeedbankClockModel();
		clockModel.initByName(
						"dormantScaling", dormantScaling,
						"tree", sbTree);
		
		// Assemble substitution model, site model
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
				
		double expResult = -41.57001439374503; // Calculated by hand
		double result = treeLikelihood.calculateLogP();
		
		System.out.println(result);
		assertEquals(expResult, result, 1e-5);
	}
	
	@Test
	public void testSeedbankClockModel3() {
		// 20 sites, dormant tips, changed branch lengths, 0.7 dormant scaling
		// Assemble data:
		Sequence sample1 = new Sequence("sample1", "TAAGAGGGGAGAGCTCTCAC");
		Sequence sample2 = new Sequence("sample2", "CTCGGAAGAGCGAGTCCCCC");
		Sequence sample3 = new Sequence("sample3", "CCCCGGTAGGACGCACCCGG");
		Sequence sample4 = new Sequence("sample4", "CCCTTATTGCCCGGGATCCA");
		
		Alignment data = new Alignment();
		data.initByName(
				"sequence", sample1,
				"sequence", sample2,
				"sequence", sample3,
				"sequence", sample4,
				"dataType", "nucleotide");
		
		// Assemble test tree:
		String newickStr =
                        "((((((sample1[&state=0]:0.3)[&state=1]:0.7,((sample2[&state=1]:0.4)"
                        + "[&state=0]:0.2)[&state=1]:0.4)[&state=1]:0.3)[&state=0]:0.3)[&state=1]:0.9,"
                        + "(sample3[&state=0]:2.4)[&state=1]:0.6)[&state=1]:1.1,sample4[&state=1]:5.0)"
                        + "[&state=1]:0.0;";
		
		SeedbankTreeFromNewick sbTreeInitializer = new SeedbankTreeFromNewick();
		SeedbankTree sbTree = new SeedbankTree();
		
		sbTreeInitializer.initByName(
                        "value", newickStr,
                        "typeLabel", "state",
                        "initial", sbTree,
                        "data", data);
		
		sbTreeInitializer.initStateNodes();
		
		// Assemble clock model
		RealParameter dormantScaling = new RealParameter();
		dormantScaling.initByName("value", "0.7");
		
		SeedbankClockModel clockModel = new SeedbankClockModel();
		clockModel.initByName(
						"dormantScaling", dormantScaling,
						"tree", sbTree);
		
		// Assemble substitution model, site model
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
				
		double expResult = -110.76258189959177; // Calculated by hand
		double result = treeLikelihood.calculateLogP();
		
		System.out.println(result);
		assertEquals(expResult, result, 1e-5);
	}
	
	@Test
	public void testSeedbankClockModel4() {
		// Single site, GTR
		// Assemble data:
		Sequence sample0 = new Sequence("sample0", "A");
		Sequence sample1 = new Sequence("sample1", "C");
		Sequence sample2 = new Sequence("sample2", "G");
		
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
		RealParameter dormantScaling = new RealParameter();
		dormantScaling.initByName("value", "0.5");
		
		SeedbankClockModel clockModel = new SeedbankClockModel();
		clockModel.initByName(
						"dormantScaling", dormantScaling,
						"tree", sbTree);
		
		// Assemble substitution model, site model
		Frequencies frequencies = new Frequencies();
		frequencies.initByName("frequencies", new RealParameter(new Double[]{1./2, 1./6, 1./4, 1./12}));
		
		GTR substModel = new GTR();
		substModel.initByName(
						"frequencies", frequencies,
						"rateAC", new Constant("0.5"),
						"rateAG", new Constant("0.75"),
						"rateAT", new Constant("1.25"),
						"rateCG", new Constant("0.75"),
						"rateCT", new Constant("1.25"),
						"rateGT", new Constant("1.5"));
	
		SiteModel siteModel = new SiteModel();
		siteModel.initByName("substModel", substModel);
		
		TreeLikelihood treeLikelihood = new TreeLikelihood();
		treeLikelihood.initByName(
						"data", data,
						"tree", sbTree,
						"siteModel", siteModel,
						"branchRateModel", clockModel,
						"implementation", "beast.evolution.likelihood.TreeLikelihood");
				
		double expResult = -3.8955794359863996; // Calculated by hand
		double result = treeLikelihood.calculateLogP();
		
		System.out.println(result);
		assertEquals(expResult, result, 1e-5);
	}
	
	@Test
	public void testSeedbankClockModel5() {
		// 10 sites
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
		RealParameter dormantScaling = new RealParameter();
		dormantScaling.initByName("value", "0.5");
		
		SeedbankClockModel clockModel = new SeedbankClockModel();
		clockModel.initByName(
						"dormantScaling", dormantScaling,
						"tree", sbTree);
		
		// Assemble substitution model, site model
		Frequencies frequencies = new Frequencies();
		frequencies.initByName("frequencies", new RealParameter(new Double[]{1./2, 1./6, 1./4, 1./12}));
		
		GTR substModel = new GTR();
		substModel.initByName(
						"frequencies", frequencies,
						"rateAC", new Constant("0.5"),
						"rateAG", new Constant("0.75"),
						"rateAT", new Constant("1.25"),
						"rateCG", new Constant("0.75"),
						"rateCT", new Constant("1.25"),
						"rateGT", new Constant("1.5"));

		SiteModel siteModel = new SiteModel();
		siteModel.initByName("substModel", substModel);
		
		TreeLikelihood treeLikelihood = new TreeLikelihood();
		treeLikelihood.initByName(
						"data", data,
						"tree", sbTree,
						"siteModel", siteModel,
						"branchRateModel", clockModel,
						"implementation", "beast.evolution.likelihood.TreeLikelihood");
				
		double expResult = -50.31314307573946; // Calculated by hand
		double result = treeLikelihood.calculateLogP();
		
		System.out.println(result);
		assertEquals(expResult, result, 1e-5);
	}
	
	@Test
	public void testSeedbankClockModel6() {
		// 20 sites, dormant tips, changed branch lengths, 0.7 dormant scaling
		// Assemble data:
		Sequence sample1 = new Sequence("sample1", "TAAGAGGGGAGAGCTCTCAC");
		Sequence sample2 = new Sequence("sample2", "CTCGGAAGAGCGAGTCCCCC");
		Sequence sample3 = new Sequence("sample3", "CCCCGGTAGGACGCACCCGG");
		Sequence sample4 = new Sequence("sample4", "CCCTTATTGCCCGGGATCCA");
		
		Alignment data = new Alignment();
		data.initByName(
				"sequence", sample1,
				"sequence", sample2,
				"sequence", sample3,
				"sequence", sample4,
				"dataType", "nucleotide");
		
		// Assemble test tree:
		String newickStr =
                        "((((((sample1[&state=0]:0.3)[&state=1]:0.7,((sample2[&state=1]:0.4)"
                        + "[&state=0]:0.2)[&state=1]:0.4)[&state=1]:0.3)[&state=0]:0.3)[&state=1]:0.9,"
                        + "(sample3[&state=0]:2.4)[&state=1]:0.6)[&state=1]:1.1,sample4[&state=1]:5.0)"
                        + "[&state=1]:0.0;";
		
		SeedbankTreeFromNewick sbTreeInitializer = new SeedbankTreeFromNewick();
		SeedbankTree sbTree = new SeedbankTree();
		
		sbTreeInitializer.initByName(
                        "value", newickStr,
                        "typeLabel", "state",
                        "initial", sbTree,
                        "data", data);
		
		sbTreeInitializer.initStateNodes();
		
		// Assemble clock model
		RealParameter dormantScaling = new RealParameter();
		dormantScaling.initByName("value", "0.7");
		
		SeedbankClockModel clockModel = new SeedbankClockModel();
		clockModel.initByName(
						"dormantScaling", dormantScaling,
						"tree", sbTree);
		
		// Assemble substitution model, site model
		Frequencies frequencies = new Frequencies();
		frequencies.initByName("frequencies", new RealParameter(new Double[]{1./2, 1./6, 1./4, 1./12}));
		
		GTR substModel = new GTR();
		substModel.initByName(
						"frequencies", frequencies,
						"rateAC", new Constant("0.5"),
						"rateAG", new Constant("0.75"),
						"rateAT", new Constant("1.25"),
						"rateCG", new Constant("0.75"),
						"rateCT", new Constant("1.25"),
						"rateGT", new Constant("1.5"));
	
		SiteModel siteModel = new SiteModel();
		siteModel.initByName("substModel", substModel);
		
		TreeLikelihood treeLikelihood = new TreeLikelihood();
		treeLikelihood.initByName(
						"data", data,
						"tree", sbTree,
						"siteModel", siteModel,
						"branchRateModel", clockModel,
						"implementation", "beast.evolution.likelihood.TreeLikelihood");
				
		double expResult = -122.77858073948788; // Calculated by hand
		double result = treeLikelihood.calculateLogP();
		
		System.out.println(result);
		assertEquals(expResult, result, 1e-5);
	}
	
	@Test
	public void testSeedbankClockModelLarge() {
		// 20 samples, 20 sites, dormant tips, changed branch lengths, 0.09 dormant scaling
		// Assemble data:
		Sequence sample0 = new Sequence("sample0", "ATAACCATGAGTTCTCGACA");
		Sequence sample1 = new Sequence("sample1", "AGGTGTGTAGCGCACCGCTT");
		Sequence sample2 = new Sequence("sample2", "GGGTACATGCCAACTAGAGA");
		Sequence sample3 = new Sequence("sample3", "ACATGCTACGATGATGTCTC");
		Sequence sample4 = new Sequence("sample4", "TGGCCCCTCTTGTTGACAGT");
		Sequence sample5 = new Sequence("sample5", "GTACTGAAACTTGTGGGTCA");
		Sequence sample6 = new Sequence("sample6", "ACACGAATGCAAGCGTAACT");
		Sequence sample7 = new Sequence("sample7", "CCTTAGCAAGTAAATGGGTC");
		Sequence sample8 = new Sequence("sample8", "CAGAGGGCAAATAACTCCGG");
		Sequence sample9 = new Sequence("sample9", "GGTTGGGGCCTAGGCAAAAG");
		Sequence sample10 = new Sequence("sample10", "AGTTAGTGTATGCGAGATTG");
		Sequence sample11 = new Sequence("sample11", "GTCCCGAGCAGACATCCTAT");
		Sequence sample12 = new Sequence("sample12", "CTCTGCCGATTCTTCACCGG");
		Sequence sample13 = new Sequence("sample13", "AGTTGTCAACTCAAGCGTGT");
		Sequence sample14 = new Sequence("sample14", "TGTGCTTTCATAACCATTTT");
		Sequence sample15 = new Sequence("sample15", "ACACAATATTTTACATGTAC");
		Sequence sample16 = new Sequence("sample16", "TCAAACCTCCGGCATCGATT");
		Sequence sample17 = new Sequence("sample17", "TGGACTCTGCAATAGCTGCG");
		Sequence sample18 = new Sequence("sample18", "GAACCACGTTCTCTGGCAAG");
		Sequence sample19 = new Sequence("sample19", "ATTATCGGTCCTATGGCAAC");
		
		Alignment data = new Alignment();
		data.initByName(
				"sequence", sample0, "sequence", sample1, "sequence", sample2, 
				"sequence", sample3, "sequence", sample4, "sequence", sample5, 
				"sequence", sample6, "sequence", sample7, "sequence", sample8, 
				"sequence", sample9, "sequence", sample10, "sequence", sample11, 
				"sequence", sample12, "sequence", sample13, "sequence", sample14, 
				"sequence", sample15, "sequence", sample16, "sequence", sample17, 
				"sequence", sample18, "sequence", sample19,
				"dataType", "nucleotide");
		
		// Assemble test tree:
		String newickStr = "((((((((((((sample0[&state=1]:1)[&state=0]:1)[&state=1]:5,"
				+ "((sample1[&state=1]:2)[&state=0]:2)[&state=1]:6)[&state=1]:1,((samp"
				+ "le2[&state=1]:3)[&state=0]:3)[&state=1]:7)[&state=1]:2,((sample3[&s"
				+ "tate=1]:4)[&state=0]:4)[&state=1]:8)[&state=1]:3,((sample4[&state=1"
				+ "]:5)[&state=0]:5)[&state=1]:9)[&state=1]:4,((sample5[&state=1]:6)[&"
				+ "state=0]:6)[&state=1]:10)[&state=1]:5,((sample6[&state=1]:7)[&state"
				+ "=0]:7)[&state=1]:11)[&state=1]:6,((sample7[&state=1]:8)[&state=0]:8"
				+ ")[&state=1]:12)[&state=1]:7,((sample8[&state=1]:9)[&state=0]:9)[&st"
				+ "ate=1]:13)[&state=1]:8,((sample9[&state=1]:10)[&state=0]:10)[&state"
				+ "=1]:14)[&state=1]:9,((((((((((sample10[&state=0]:10)[&state=1]:10,("
				+ "sample11[&state=0]:11)[&state=1]:11)[&state=1]:1,(sample12[&state=0"
				+ "]:12)[&state=1]:12)[&state=1]:2,(sample13[&state=0]:13)[&state=1]:1"
				+ "3)[&state=1]:3,(sample14[&state=0]:14)[&state=1]:14)[&state=1]:4,(s"
				+ "ample15[&state=0]:15)[&state=1]:15)[&state=1]:5,(sample16[&state=0]"
				+ ":16)[&state=1]:16)[&state=1]:6,(sample17[&state=0]:17)[&state=1]:17"
				+ ")[&state=1]:7,(sample18[&state=0]:18)[&state=1]:18)[&state=1]:8,(sa"
				+ "mple19[&state=0]:19)[&state=1]:19)[&state=1]:9)[&state=1]:0;";
		
		SeedbankTreeFromNewick sbTreeInitializer = new SeedbankTreeFromNewick();
		SeedbankTree sbTree = new SeedbankTree();
		
		sbTreeInitializer.initByName(
                        "value", newickStr,
                        "typeLabel", "state",
                        "initial", sbTree,
                        "data", data);
		
		sbTreeInitializer.initStateNodes();
		
		// Assemble clock model
		RealParameter dormantScaling = new RealParameter();
		dormantScaling.initByName("value", "0.09");
		
		SeedbankClockModel clockModel = new SeedbankClockModel();
		clockModel.initByName(
						"dormantScaling", dormantScaling,
						"tree", sbTree);
		
		// Assemble substitution model
		Frequencies frequencies = new Frequencies();
		frequencies.initByName("frequencies", new RealParameter(new Double[]{1./2, 1./6, 1./4, 1./12}));
		
		GTR substModel = new GTR();
		substModel.initByName(
						"frequencies", frequencies,
						"rateAC", new Constant("0.5"),
						"rateAG", new Constant("0.75"),
						"rateAT", new Constant("1.25"),
						"rateCG", new Constant("0.75"),
						"rateCT", new Constant("1.25"),
						"rateGT", new Constant("1.5"));
	
		SiteModel siteModel = new SiteModel();
		siteModel.initByName("substModel", substModel);
		
		TreeLikelihood treeLikelihood = new TreeLikelihood();
		treeLikelihood.initByName(
						"data", data,
						"tree", sbTree,
						"siteModel", siteModel,
						"branchRateModel", clockModel,
						"implementation", "beast.evolution.likelihood.TreeLikelihood");
				
		double expResult = -637.4025254166279; // Calculated by hand
		double result = treeLikelihood.calculateLogP();
		
		System.out.println(result);
		assertEquals(expResult, result, 1e-5);
	}
}

/////
//clockModel.calculateRates(sbTree.getRoot());
//
//for (Node n : sbTree.getNodesAsArray()) {
//	SeedbankNode sbNode = (SeedbankNode) n;
//	System.out.println(String.format("ID: %s, type: %d, height: %f", sbNode.getID(), sbNode.getNodeType(), sbNode.getHeight()));
//	System.out.println(clockModel.getRateForBranch(sbNode));
//}
/////
