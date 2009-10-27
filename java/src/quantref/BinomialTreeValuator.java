package quantref;

import org.joda.time.Days;

/**
 * User: derekbennett
 * Date: Oct 24, 2009
 * Time: 3:39:49 PM
 */
public class BinomialTreeValuator implements Valuator<Option>{
    private StateOfTheWorld sotw;
    private int N=200;

    public BinomialTreeValuator(StateOfTheWorld sotw, int numberOfTreeSteps) {
        this.sotw = sotw;
        this.N = numberOfTreeSteps;
    }

    public void setStateOfTheWorld(StateOfTheWorld sotw) {
        this.sotw =sotw;
    }

    private double binomialTreeValue(Poc poc, OptionType type, double So, double sigma, double tte, double K, double r, int N){
        double dt   = tte/N;
        double nu   = r - 0.5*Math.pow(sigma,2);
        double dxu  = Math.sqrt(Math.pow(sigma,2)*dt + Math.pow(nu*dt,2));
        double dxd  = -dxu;
        double pu   = 0.5 + 0.5*(nu*dt/dxu);
        double pd   = 1 - pu;
        double disc = Math.exp(-r*dt);
        double dpu = disc*pu;
        double dpd = disc*pd;
        double edxud = Math.exp(dxu-dxd);
        double edxd = Math.exp(dxd);

        //Initialize asset prices at maturity
        double [] St = new double[N+1];
        St[0] = So*Math.exp(N*dxd);
        for (int i=1;i<=N;i++) {
            St[i] = St[i-1]*edxud;
        }

        //Initialize option values at maturity
        double [] V = new double[N+1];
        for (int j=0;j<=N;j++) {
            V[j] = poc.equals(Poc.CALL) ? Math.max(0.0,St[j]-K) : Math.max(0.0,K-St[j]);
        }

        //Step backwards in time to now to get the initial option price
        for (int i=N; i>=0; i--) {
            for(int j=0;j<i;j++) {
                V[j] = dpu*V[j+1] + dpd*V[j];
                if (type.equals(OptionType.AMERICAN)) {
                    St[j] = St[j]/edxd;
                    V[j] = poc.equals(Poc.CALL) ? Math.max(V[j],St[j]-K) : Math.max(V[j],K-St[j]);
                }
            }
        }
        return(V[0]);
    }

    public Valuation getValuation(Option option) {
        double So = sotw.getPrices().get(option.getUnderlying());
        double K = option.getStrike();
        double r = sotw.getInterestRate();
        double sigma = sotw.getVolatilities().get(option);
        double tte = Days.daysBetween(sotw.getBusinessDate(),option.getExpiry()).getDays() / 365;
        double price = 0;
        price = binomialTreeValue(option.getPoc(),option.getType(),So,sigma,tte,K,r,N);
        return new Valuation<Option>(option, price);

    }

    public Valuation getImpliedVolatility(Option instrument, double price) {
        return null;
    }
}
