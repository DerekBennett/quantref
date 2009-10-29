package quantref;

import org.joda.time.Days;

/**
 * User: derekbennett
 * Date: Oct 24, 2009
 * Time: 3:39:49 PM
 */
public class TrinomialTreeValuator implements Valuator<Option>{
    private StateOfTheWorld sotw;
    private int N=200;

    public TrinomialTreeValuator(StateOfTheWorld sotw, int numberOfTreeSteps) {
        this.sotw = sotw;
        this.N = numberOfTreeSteps;
    }

    public void setStateOfTheWorld(StateOfTheWorld sotw) {
        this.sotw =sotw;
    }

    private double trinomialTreeValue(Poc poc, OptionType type, double So, double sigma, double tte,
                                      double K, double r, double div, int N){
        double dt   = tte/N;
        double dx   = sigma*Math.sqrt(3*dt);
        double nu   = r - div - 0.5*Math.pow(sigma,2);
        double edx  = Math.exp(dx);
        double xi   = (Math.pow(sigma,2)*dt + Math.pow(nu,2)*Math.pow(dt,2))/Math.pow(dx,2);
        double pu   = (xi + nu*dt/dx)/2;
        double pm   = 1.0 - xi;
        double pd   = (xi - nu*dt/dx)/2;
        double disc = Math.exp(-r*dt);

        //Initialize asset prices at maturity
        double [] St = new double[2*N+1];
        St[0] = So*Math.exp(-N*dx);
        for (int j=1;j<=2*N;j++) {
            St[j] = St[j-1]*edx;
        }

        //Initialize option values at maturity
        double [][] V = new double[N+1][2*N+1];
        for (int j=0;j<=2*N;j++) {
            V[N][j] = poc.equals(Poc.CALL) ? Math.max(0.0,St[j]-K) : Math.max(0.0,K-St[j]);
        }

        //Step backwards in time to now to get the initial option price
        int treeNarrowsBy=0;
        for (int i=N-1; i>=0; i--) {
            treeNarrowsBy++;
            for(int j=treeNarrowsBy;j<=2*N-treeNarrowsBy;j++) {
                V[i][j] = disc*(pu*V[i+1][j+1] + pm*V[i+1][j] + pd*V[i+1][j-1]);
                if (type.equals(OptionType.AMERICAN)) {
                    St[j] = St[j]/edx;//This adjustment is wrong
                    V[i][j] = poc.equals(Poc.CALL) ? Math.max(V[i][j],St[j]-K) : Math.max(V[i][j],K-St[j]);
                }
            }
        }
        return(V[0][N]);
    }

    public Valuation getValuation(Option option) {
        double So = sotw.getPrices().get(option.getUnderlying());
        double K = option.getStrike();
        double r = sotw.getInterestRate();
        double div = option.getUnderlying().getContinuousDividend();
        double sigma = sotw.getVolatilities().get(option);
        double tte = Days.daysBetween(sotw.getBusinessDate(),option.getExpiry()).getDays() / 365;
        double price = 0;
        price = trinomialTreeValue(option.getPoc(),option.getType(),So,sigma,tte,K,r,div,N);
        return new Valuation<Option>(option, price);

    }

    public Valuation getImpliedVolatility(Option instrument, double price) {
        return null;
    }
}