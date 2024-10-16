package seedbanktree.inference.distribution;

import org.apache.commons.math.MathException;
import org.apache.commons.math.MathRuntimeException;
import org.apache.commons.math.distribution.Distribution;
import org.apache.commons.math.distribution.IntegerDistribution;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.distribution.ParametricDistribution;

public class Bernoulli extends ParametricDistribution {
	final public Input<Function> pInput = new Input<>("p", "bernoulli success probability", Validate.REQUIRED);
	BernoulliImpl dist = new BernoulliImpl();

    @Override
    public void initAndValidate() {
    	refresh();
    }
    
    /**
     * ensure internal state is up to date *
     */
    void refresh() {
        double p = pInput.get().getArrayValue();
        dist.setP(p);
    }

    @Override
    public Distribution getDistribution() {
    	refresh();
        return dist;
    }

    class BernoulliImpl implements IntegerDistribution {
    	double p;
    	
    	public void setP(double p) {
    		this.p = p;
    	}
    	
		@Override
		public double probability(double x) {
			if (x == 1) return p;
			else if (x == 0) return 1-p;
			else return 0;
		}

		@Override
		public double cumulativeProbability(double x) throws MathException {
			if (x >= 1) return 1;
			else if (x >= 0) return 1-p;
			else return 0;
		}

		@Override
		public double cumulativeProbability(double x0, double x1) throws MathException {
			if (x0 > x1) {
	            throw MathRuntimeException.createIllegalArgumentException(
                    "lower endpoint ({0}) must be less than or equal to upper endpoint ({1})",
                    x0, x1);
	        }
	        return cumulativeProbability(x1) - cumulativeProbability(x0);
		}

		@Override
		public double probability(int x) {
			return probability((double) x);
		}

		@Override
		public double cumulativeProbability(int x) throws MathException {
			return cumulativeProbability((double) x);
		}

		@Override
		public double cumulativeProbability(int x0, int x1) throws MathException {
			return cumulativeProbability((double) x0, (double) x1);
		}

		@Override
		public int inverseCumulativeProbability(double p) throws MathException {
			if (this.p == 1) {
				if (p < 1) return 0;
				else return Integer.MAX_VALUE;
			} else if (this.p > 0) {				
				if (p < (1-this.p)) return -1;
				else if (p < 1) return 0;
				else return Integer.MAX_VALUE;
			} else { // this.p == 0
				if (p < 1) return -1;
				else return Integer.MAX_VALUE;
			}
		}

    } 

}
