package seedbanktree.util;

import java.io.PrintStream;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;
import beast.base.inference.CalculationNode;

public class CKLogger extends CalculationNode implements Loggable {

	final public Input<Function> cInput = new Input<>("c", "c, transition rate.", Validate.REQUIRED);
	final public Input<Function> KInput = new Input<>("K", "K, relative seedbank population size.", Validate.REQUIRED);
	
	@Override
	public void initAndValidate() { }
	
	@Override
	public void init(PrintStream out) {
        out.print("c*K\t");
	}

	@Override
	public void log(long sample, PrintStream out) {
        out.print(String.valueOf(cInput.get().getArrayValue() * KInput.get().getArrayValue()) + "\t");
	}

	@Override
	public void close(PrintStream out) {

	}
}
