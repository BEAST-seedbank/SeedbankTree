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
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;
import cern.colt.Arrays;

public class LambdaIndicatorLogger extends CalculationNode implements Loggable, Function {
	// i.e. lambdas
	final public Input<Function> lambdasInput = new Input<>("lambdas", "", Validate.REQUIRED);
	final public Input<Function> indicatorsInput = new Input<>("indicators", "", Validate.REQUIRED);
//	final public Input<RandomLocalClockModel> RLCInput = new Input<>("RLC", "", Validate.REQUIRED);
	final public Input<Tree> treeInput = new Input<>("tree", "");
	
	int[] indicatorsSum;
	double[] lambdasSum;
	
	double[] indicatorsAvg;
	double[] lambdasAvg;
	
	double count = 0;

	@Override
	public void initAndValidate() {
	}

	public void update() {
	}

	@Override
	public void init(PrintStream out) {
		out.print("indicatorsAvg\t");
		out.print("lambdasAvg\t");
		
		int dim = lambdasInput.get().getDimension();
		indicatorsSum = new int[dim];
		lambdasSum = new double[dim];
		indicatorsAvg = new double[dim];
		lambdasAvg = new double[dim];
		
	}

	@Override
	public void log(long sample, PrintStream out) {
		count += 1;
		
		Function indicators = indicatorsInput.get();
		Function lambdas = lambdasInput.get();
		int numActive = 0;
		for (int i = 0; i < lambdas.getDimension(); i++) {
			if (indicators.getArrayValue(i) == 1) {
				assert (lambdas.getArrayValue(i) != 0.0);
				numActive++;
				
				indicatorsSum[i] += 1;
				lambdasSum[i] += lambdas.getArrayValue(i);
				indicatorsAvg[i] = indicatorsSum[i] / count;
				lambdasAvg[i] = lambdasSum[i] / count;
			}
		}

		out.print(numActive + "\t");
		out.print(Arrays.toString(indicatorsAvg) + "\t");
		out.print(Arrays.toString(lambdasAvg) + "\t");
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
