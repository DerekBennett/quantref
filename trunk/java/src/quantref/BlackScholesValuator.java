package quantref;

import org.joda.time.Days;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.MathException;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.ConvergenceException;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolver;
import org.apache.commons.math.analysis.solvers.UnivariateRealSolverFactory;
import org.apache.commons.math.analysis.UnivariateRealFunction;

/**
 * User: derekbennett
 * Date: Oct 18, 2009
 * Time: 10:14:18 PM
 */
public class BlackScholesValuator implements Valuator<Option> {
    private StateOfTheWorld sotw;

    public void setStateOfTheWorld(final StateOfTheWorld sotw) {
        this.sotw =sotw;
    }

    public Valuation getValuation(final Option option) {
        double So = sotw.getPrices().get(option.getUnderlying());
        double K = option.getStrike();
        double r = sotw.getInterestRate();
        double sigma = sotw.getVolatilities().get(option);
        double tte = Days.daysBetween(sotw.getBusinessDate(),option.getExpiry()).getDays() / 365.0;
        double price = 0;
        try {
            price = blackScholesValue(option.getPoc(),So,sigma,tte,K,r);
            return new Valuation<Option>(option, price);
        } catch (MathException e) {
            e.printStackTrace();
            return null;
        }
    }

    private double blackScholesValue(Poc poc, double So, double sigma, double tte, double K, double r)
            throws MathException{
        double logPriceRatio = Math.log(So / K);
        double rateTimeVol = r + (0.5 * Math.pow(sigma, 2) * tte);
        double sigmaRootTte = sigma * Math.sqrt(tte);
        double d1 = (logPriceRatio + rateTimeVol) / sigmaRootTte;
        double d2 = d1 - sigmaRootTte;
        NormalDistributionImpl N = new NormalDistributionImpl();
        double price = 0;
        if (poc.equals(Poc.PUT)) {
            price = K * Math.exp(-r * tte) * N.cumulativeProbability(-d2) - So * N.cumulativeProbability(-d1);
        } else {
            price = So * N.cumulativeProbability(d1) - K * Math.exp(-r * tte) * N.cumulativeProbability(d2);
        }
        return price;
    }

    public Valuation getImpliedVolatility(final Option option, final double price) {
        final double So = sotw.getPrices().get(option.getUnderlying());
        final double K = option.getStrike();
        final double r = sotw.getInterestRate();
        final double tte = Days.daysBetween(sotw.getBusinessDate(),option.getExpiry()).getDays() / 365;
        UnivariateRealSolverFactory factory = UnivariateRealSolverFactory.newInstance();
        UnivariateRealSolver rootFinder = factory.newBisectionSolver();
        UnivariateRealFunction blackScholes = new UnivariateRealFunction() {
            public double value(double impliedVol) throws FunctionEvaluationException {
                try {
                    return blackScholesValue(option.getPoc(),So,impliedVol,tte,K,r)-price;
                } catch (MathException e) {
                    e.printStackTrace();
                    return 0;
                }
            }
        };
        try {
            return new Valuation<Option>(option,rootFinder.solve(blackScholes,0,1));
        } catch (ConvergenceException e) {
            e.printStackTrace();
        } catch (FunctionEvaluationException e) {
            e.printStackTrace();
        }
        return new Valuation<Option>(option,0.0);
    }
}
