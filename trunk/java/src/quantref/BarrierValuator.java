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
 * WARNING: This class currently only implements an American up-and-in Call option
 *          The plan is to add some kind of pay-off lambda function that can be passed into these valuators 
 */
public class BarrierValuator implements Valuator<Option>{
    private StateOfTheWorld sotw;
    private int N=200;
    private int M=1000;
    public double stddev=0.0;
    public double stderr=0.0;
    private double H=100;

    public BarrierValuator(StateOfTheWorld sotw, int numTimeSteps, int numPaths, double barrierPrice) {
        this.sotw = sotw;
        this.N = numTimeSteps;
        this.M = numPaths;
        this.H = barrierPrice;
    }

    public void setStateOfTheWorld(StateOfTheWorld sotw) {
        this.sotw =sotw;
    }

    public Valuation getValuation(Option option) {
        Poc        poc  = option.getPoc();
        OptionType type = option.getType();
        double So = sotw.getPrices().get(option.getUnderlying());
        double K = option.getStrike();
        double r = sotw.getInterestRate();
        double div = option.getUnderlying().getContinuousDividend();
        double sigma = sotw.getVolatilities().get(option);
        double tte = Days.daysBetween(sotw.getBusinessDate(),option.getExpiry()).getDays() / 365.0;
        double dt   = tte/N;
        double nu   = r - div - 0.5*Math.pow(sigma,2);
        double nudt = nu*dt;
        double sigsdt= sigma*Math.sqrt(dt);

        Random random = new Random();
        double sum_Vt =0.0;
        double sum_Vt2=0.0;
        double St;
        for (int j=0;j<M;j++) {
            // For each iteration, we want to keep track of option prices at each step
            // if the barrier is crossed, get the optimal exericse after that crossing
            St = So;
            double [] Vt = new double[N];
            double actualVt = 0.0;
            int crossedBarrierAt=0;
            int i;
            for (i=0;i<N;i++) {
                St = St *Math.exp(nudt + sigsdt*random.nextGaussian());
                if (St>H){
                    crossedBarrierAt=i;
                }
                Vt[i] = poc.equals(Poc.CALL) ? Math.max(0.0,St-K) : Math.max(0.0,K-St);
            }
            if (crossedBarrierAt>0){
                for (int x=crossedBarrierAt; x<N; x++){
                    actualVt = Math.max(actualVt,Vt[x]);
                }
            }
            else {
                actualVt=0;
            }

            sum_Vt  = sum_Vt + actualVt;
            sum_Vt2 = sum_Vt2 + actualVt*actualVt;
        }
        this.stddev = Math.sqrt((sum_Vt2 - sum_Vt*sum_Vt/M)*Math.exp(-2*r*tte)/M-1);
        this.stderr = stddev/Math.sqrt(M);
        return new Valuation<Option>(option, sum_Vt/M*Math.exp(-r*tte));

    }

    public Valuation getImpliedVolatility(Option instrument, double price) {
        return null;
    }
}