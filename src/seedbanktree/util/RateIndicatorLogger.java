package seedbanktree.util;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import beast.base.evolution.branchratemodel.RandomLocalClockModel;
import beast.base.evolution.branchratemodel.RandomLocalClockModel2;
import beast.base.inference.CalculationNode;

public class RateIndicatorLogger extends CalculationNode implements Loggable, Function {
	final public Input<Function> ratesInput = new Input<>("rates", "", Validate.REQUIRED);
	final public Input<Function> indicatorsInput = new Input<>("indicators", "", Validate.REQUIRED);
	final public Input<RandomLocalClockModel> RLCInput = new Input<>("RLC", "", Validate.REQUIRED);

	// counts[0] = dormant; counts[1] = active
//    int[] counts;

	@Override
	public void initAndValidate() {
//		counts = new int[]{0, 0};
	}

	public void update() {
//		final SeedbankTree sbTree = seedbankTreeInput.get();
//		counts[0] = 0;
//		counts[1] = 0;
//		
//        for (Node node : sbTree.getNodesAsArray()) {
//        	SeedbankNode sbNode = (SeedbankNode)node;
//        	counts[sbNode.getNodeType()] += 1;
//        			
//            if (node.isRoot()) {
//                continue;
//            }
//            
//            if (sbNode.getNodeType() == 0) {
////            	assert(sbNode.getChangeCount() % 2 == 1);
//            	counts[0] += sbNode.getChangeCount() / 2;
//    			counts[1] += 1 + sbNode.getChangeCount() / 2;
//            } else {
////            	assert(sbNode.getChangeCount() % 2 == 0);
//    			counts[0] += sbNode.getChangeCount() / 2;
//    			counts[1] += sbNode.getChangeCount() / 2;
//            }
//        }
	}

	@Override
	public void init(PrintStream out) {
		out.print("rates\t");
		out.print("indicators\t");
	}

	@Override
	public void log(long sample, PrintStream out) {
		Function rates = ratesInput.get();
		Function indicators = indicatorsInput.get();
		Boolean activeRate = false;
		for (int i = 0; i < rates.getDimension(); i++) {
			if (rates.getArrayValue(i) != 1.0 && indicators.getArrayValue(i) == 0) {
				activeRate = true;
			}
		}

		out.print(activeRate + "\t");
		out.print(rates + "\t");
		out.print(indicators + "\t");

		FileWriter fw;
		try {
			fw = new FileWriter("singledormancy/manualLogging.txt", true);
			BufferedWriter bw = new BufferedWriter(fw);
			PrintWriter writer = new PrintWriter(bw);
			writer.println(activeRate + "\t");
			writer.println(rates + "\t");
			writer.println(indicators + "\t");
			writer.println("\n");
			double[] normalizedRates = RLCInput.get().getRates();
			writer.println(normalizedRates.length);
			for (int i = 0; i < normalizedRates.length; i++) {
				writer.print(normalizedRates[i] + " ");
			}
			writer.println("\n\n");
			writer.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	@Override
	public void close(PrintStream out) {

	}

	// function implementation

	@Override
	public int getDimension() {
		return 2;
	}

	@Override
	public double getArrayValue(int dim) {
		update();

		int x = 1 / 0;
		return 1.0;
	}

}
