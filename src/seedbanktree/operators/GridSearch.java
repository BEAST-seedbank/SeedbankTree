package seedbanktree.operators;

import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.util.InputUtil;

public class GridSearch extends Operator {
	
	public final Input<RealParameter> cInput = new Input<>("c", "c to search", 
			Input.Validate.REQUIRED);
	public final Input<Double> cLowerInput = new Input<>("cLower", "Lower limit for c");
	public final Input<Double> cUpperInput = new Input<>("cUpper", "Upper limit for c");
	
	public final Input<RealParameter> KInput = new Input<>("K", "K to search", 
			Input.Validate.REQUIRED);
	public final Input<Double> KLowerInput = new Input<>("KLower", "Lower limit for K");
	public final Input<Double> KUpperInput = new Input<>("KUpper", "Upper limit for K");
	
	public final Input<Integer> quantityInput = new Input<>("quantity", "Number of values to search (min 2)");
	
	Double cLower, cUpper, KLower, KUpper;
	Integer quantity;
	Double cStep, KStep;
	Double cNext, KNext;
	
	@Override
	public void initAndValidate() {
		cLower = cLowerInput.get();
		cUpper = cUpperInput.get();
		KLower = KLowerInput.get();
		KUpper = KUpperInput.get();
		quantity = quantityInput.get();
		
		cStep = (cUpper - cLower) / (quantity - 1);
		KStep = (KUpper - KLower) / (quantity - 1);
		
		cNext = cLower;
		KNext = KLower;
	}
	
	@Override
	public double proposal() {
		if (cNext > cUpper + 1e-5) return 0.0;
		if (KNext > KUpper + 1e-5) KNext = KLower;
		
		final RealParameter c = (RealParameter)InputUtil.get(cInput, this);
		final RealParameter K = (RealParameter)InputUtil.get(KInput, this);
		
		c.setValue(cNext);
		K.setValue(KNext);
		
		KNext += KStep;
		if (KNext > KUpper + 1e-5) cNext += cStep;
//		System.out.println(cNext);
//		System.out.println(KNext);
		
		return Double.POSITIVE_INFINITY;
	}

}
