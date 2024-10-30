package seedbanktree.inference.distribution;

import org.apache.commons.math.MathException;
import org.apache.commons.math.MathRuntimeException;
import org.apache.commons.math.distribution.ContinuousDistribution;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.distribution.Prior;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;


public class SpikeAndSlabPrior extends Prior {
	
	final public Input<Function> spikeLocInput = new Input<>("spike", "spike location, defaults to 0");
    final public Input<ContinuousDistribution> slabDistInput = 
    		new Input<>("slab", "the distribution that comprises the slab", Validate.REQUIRED); 
	final public Input<Function> etasInput = new Input<>("etas", "mixture weight", Validate.REQUIRED);
	
	protected SpikeAndSlabImpl spikeAndSlab = new SpikeAndSlabImpl(0.0);
	
	public SpikeAndSlabPrior() {
		distInput.setRule(Validate.FORBIDDEN);
	}
	
	@Override
    public void initAndValidate() {
		if (spikeLocInput.get() != null) {
			spikeAndSlab.setSpike(spikeLocInput.get().getArrayValue());
		}
		spikeAndSlab.setSlab(slabDistInput.get());
		
        if (m_x.get().getDimension() != etasInput.get().getDimension()) {
            throw new IllegalArgumentException(String.format("Dimension of m_x (%d) and gamma (%d) must match.", m_x.get().getDimension(), etasInput.get().getDimension()));
        }
        
        calculateLogP();
    }
	
	@Override
    public double calculateLogP() {
        Function x = m_x.get();
        Function gamma = etasInput.get();
        if (x instanceof RealParameter || x instanceof IntegerParameter) {
            // test that parameter is inside its bounds
            double l = 0.0;
            double h = 0.0;
            if (x instanceof RealParameter) {
                l = ((RealParameter) x).getLower();
                h = ((RealParameter) x).getUpper();
            } else {
                l = ((IntegerParameter) x).getLower();
                h = ((IntegerParameter) x).getUpper();
            }
            
            for (int i = 0; i < x.getDimension(); i++) {
                double value = x.getArrayValue(i);
                if (value < l || value > h) {
                    logP = Double.NEGATIVE_INFINITY;
                    System.out.println("OUT OF BOUNDS " + value);
                    return Double.NEGATIVE_INFINITY;
                }
            }
        }
        logP = 0;
        for (int i = 0; i < x.getDimension(); i++) {
            final double x_val = x.getArrayValue(i);
            final int gamma_val = (int) gamma.getArrayValue(i);
            logP += spikeAndSlab.logDensity(x_val, gamma_val);
        }
        if (logP == Double.POSITIVE_INFINITY) {
            logP = Double.NEGATIVE_INFINITY;
        }
        
        return logP;
    }
	
	class SpikeAndSlabImpl {
		double spike;
		ContinuousDistribution slab;
		
		SpikeAndSlabImpl(double spike) {
            setSpike(spike);
        }
        
        public void setSpike(double spike) {
            this.spike = spike;
        }
        
        public void setSlab(ContinuousDistribution slab) {
        	this.slab = slab;
        }
        
		public double cumulativeProbability(double x, int eta) throws MathException {
			if (eta!= 0 && eta != 1 )
				throw MathRuntimeException.createIllegalArgumentException("eta can only be 0 or 1");
			
			if (eta == 1) {
				return slab.cumulativeProbability(x);
			} else {
				return x >= spike ? 1 : 0;
			}
		}

		public double cumulativeProbability(double x0, double x1, int eta) throws MathException {
	        if (x0 > x1) {
	            throw MathRuntimeException.createIllegalArgumentException(
                    "lower endpoint ({0}) must be less than or equal to upper endpoint ({1})",
                    x0, x1);
	        }
	        return cumulativeProbability(x1, eta) - cumulativeProbability(x0, eta);
	    }

		public double density(double x, int eta) {
			if (eta!= 0 && eta != 1 )
				throw MathRuntimeException.createIllegalArgumentException("eta can only be 0 or 1");
			
			assert eta == 0 || eta == 1;
			
			if (eta == 1) {
				return slab.density(x);
			} else {
				return x == spike ? 1 : 0;
			}
		}

		public double logDensity(double x, int eta) {
			if (eta!= 0 && eta != 1 )
				throw MathRuntimeException.createIllegalArgumentException("eta can only be 0 or 1");
			
			if (eta == 1) {
				return slab.logDensity(x);
			} else {
				return x == spike ? 0 : Double.NEGATIVE_INFINITY;
			}
		}
	}
}