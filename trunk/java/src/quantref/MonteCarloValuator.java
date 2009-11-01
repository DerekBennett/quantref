package quantref;

import org.joda.time.Days;
import org.apache.commons.math.random.RandomDataImpl;
import org.apache.commons.math.random.RandomData;
import org.apache.commons.math.random.JDKRandomGenerator;

import java.util.Random;

/**
 * User: derekbennett
 * Date: Oct 24, 2009
 * Time: 3:39:49 PM
 */
public class MonteCarloValuator implements Valuator<Option>{
    private StateOfTheWorld sotw;
    private int N=200;
    private int M=1000;
    public double stddev=0.0;
    public double stderr=0.0;

    public MonteCarloValuator(StateOfTheWorld sotw, int numTimeSteps, int numPaths) {
        this.sotw = sotw;
        this.N = numTimeSteps;
        this.M = numPaths;
    }

    public void setStateOfTheWorld(StateOfTheWorld sotw) {
        this.sotw =sotw;
    }

    private double monteCarloValue(Poc poc, OptionType type, double So, double sigma, double tte,
                                      double K, double r, double div, int N, int M){
        double dt   = tte/N;
        double nu   = r - div - 0.5*Math.pow(sigma,2);
        double nudt = nu*dt;
        double sigsdt= sigma*Math.sqrt(dt);
        double lnS   = Math.log(So);

        Random random = new Random();
        double lnSt;
        double sum_Vt =0.0;
        double sum_Vt2=0.0;
        double St;
        double Vt;
        for (int j=0;j<M;j++) {
            lnSt = lnS;
            for (int i=0;i<N;i++) {
                lnSt = lnSt + nudt + sigsdt*random.nextGaussian();
            }
            St = Math.exp(lnSt);
            Vt = poc.equals(Poc.CALL) ? Math.max(0.0,St-K) : Math.max(0.0,K-St);
            sum_Vt  = sum_Vt + Vt;
            sum_Vt2 = sum_Vt2 + Vt*Vt;
        }
        this.stddev = Math.sqrt((sum_Vt2 - sum_Vt*sum_Vt/M)*Math.exp(-2*r*tte)/M-1);
        this.stderr = stddev/Math.sqrt(M);
        return sum_Vt/M*Math.exp(-r*tte);
    }

    public Valuation getValuation(Option option) {
        double So = sotw.getPrices().get(option.getUnderlying());
        double K = option.getStrike();
        double r = sotw.getInterestRate();
        double div = option.getUnderlying().getContinuousDividend();
        double sigma = sotw.getVolatilities().get(option);
        double tte = Days.daysBetween(sotw.getBusinessDate(),option.getExpiry()).getDays() / 365.0;
        double price = 0;
        price = monteCarloValue(option.getPoc(),option.getType(),So,sigma,tte,K,r,div,N, M);
        return new Valuation<Option>(option, price);

    }

    public Valuation getImpliedVolatility(Option instrument, double price) {
        return null;
    }
}