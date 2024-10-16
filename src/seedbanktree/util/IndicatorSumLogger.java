package seedbanktree.util;

import java.io.PrintStream;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Tree;
import beast.base.inference.CalculationNode;

public class IndicatorSumLogger extends CalculationNode implements Loggable, Function {
	final public Input<Function> indicatorsInput = new Input<>("indicators", "", Validate.REQUIRED);
	final public Input<Tree> treeInput = new Input<>("tree", "");
	
	int sum;
	
	@Override
	public void initAndValidate() {}

	@Override
	public void init(PrintStream out) {
		out.print("indicatorsSum\t");
	}

	@Override
	public void log(long sample, PrintStream out) {
		Function indicators = indicatorsInput.get();
		sum = 0;
		for (int i = 0; i < indicators.getDimension(); i++)
			sum += indicators.getArrayValue(i);

		out.print(sum + "\t");
	}

	@Override
	public void close(PrintStream out) {

	}

	// Function implementation

	@Override
	public int getDimension() {
		return 1;
	}

	@Override
	public double getArrayValue(int dim) {
		return sum;
	}

}
