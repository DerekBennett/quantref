package quantref;

import org.joda.time.Days;

/**
 * User: derekbennett
 * Date: Oct 24, 2009
 * Time: 3:39:49 PM
 */
public class ExplicitFiniteDifferenceValuator implements Valuator<Option>{
    private StateOfTheWorld sotw;
    private int N=200;
    private int Nj=200;

    public ExplicitFiniteDifferenceValuator(StateOfTheWorld sotw, int timeDimension, int spaceDimension) {
        this.sotw = sotw;
        this.N  = timeDimension;
        this.Nj = spaceDimension;
    }

    public void setStateOfTheWorld(StateOfTheWorld sotw) {
        this.sotw =sotw;
    }

    private double explicitFiniteDifferenceValue(Poc poc, OptionType type, double So, double sigma, double tte,
                                      double K, double r, double div, int N, int Nj){
        double dt   = tte/N;
        double dx   = sigma*Math.sqrt(3*dt);
        double nu   = r - div - 0.5*Math.pow(sigma,2);
        double edx  = Math.exp(dx);
        double xi   = Math.pow(sigma/dx,2);
        double pu   = 0.5*dt*(xi + nu/dx);
        double pm   = 1.0 - dt*xi - r*dt;
        double pd   = 0.5*dt*(xi - nu/dx);

        //Initialize asset prices at maturity
        double [] St = new double[2*Nj+1];
        St[0] = So*Math.exp(-N*dx);
        for (int j=1;j<=2*Nj;j++) {
            St[j] = St[j-1]*edx;
        }

        //Initialize option values at maturity
        double [][] V = new double[N+1][2*Nj+1];
        for (int j=0;j<=2*Nj;j++) {
            V[N][j] = poc.equals(Poc.CALL) ? Math.max(0.0,St[j]-K) : Math.max(0.0,K-St[j]);
        }

        //Step backwards in time to now to get the initial option price
        for (int i=N-1; i>=0; i--) {
            for(int j=1;j<=2*Nj-1;j++) {
                V[i][j] = pu*V[i+1][j+1] + pm*V[i+1][j] + pd*V[i+1][j-1]; // TODO doesn't discount here (unlike trinomial)
                if (type.equals(OptionType.AMERICAN)) {
                    V[i][j] = poc.equals(Poc.CALL) ? Math.max(V[i][j],St[j]-K) : Math.max(V[i][j],K-St[j]);
                }
            }
            V[i][0]=V[i][1];
            V[i][2*Nj]=V[i][2*Nj-1] + (St[2*Nj]-St[2*Nj-1]);
        }
        return(V[0][Nj]);
    }

    public Valuation getValuation(Option option) {
        double So = sotw.getPrices().get(option.getUnderlying());
        double K = option.getStrike();
        double r = sotw.getInterestRate();
        double div = option.getUnderlying().getContinuousDividend();
        double sigma = sotw.getVolatilities().get(option);
        double tte = Days.daysBetween(sotw.getBusinessDate(),option.getExpiry()).getDays() / 365.0;
        double price = 0;
        price = explicitFiniteDifferenceValue(option.getPoc(),option.getType(),So,sigma,tte,K,r,div,N,Nj);
        return new Valuation<Option>(option, price);

    }

    public Valuation getImpliedVolatility(Option instrument, double price) {
        return null;
    }
}