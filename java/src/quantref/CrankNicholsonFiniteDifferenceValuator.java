package quantref;

import org.joda.time.Days;

/**
 * User: derekbennett
 * Date: Nov 1, 2009
 * Time: 7:05:49 PM
 */
public class CrankNicholsonFiniteDifferenceValuator implements Valuator<Option>{
    private StateOfTheWorld sotw;
    private int N=200;
    private int Nj=200;

    public CrankNicholsonFiniteDifferenceValuator(StateOfTheWorld sotw, int timeDimension, int spaceDimension) {
        this.sotw = sotw;
        this.N  = timeDimension;
        this.Nj = spaceDimension;
    }

    public void setStateOfTheWorld(StateOfTheWorld sotw) {
        this.sotw =sotw;
    }

    public Valuation getValuation(Option option) {
        Poc poc         = option.getPoc();
        OptionType type = option.getType();
        double So   = sotw.getPrices().get(option.getUnderlying());
        double K    = option.getStrike();
        double r    = sotw.getInterestRate();
        double div  = option.getUnderlying().getContinuousDividend();
        double sigma= sotw.getVolatilities().get(option);
        double tte  = Days.daysBetween(sotw.getBusinessDate(),option.getExpiry()).getDays() / 365.0;
        double dt   = tte/N;
        double dx   = sigma*Math.sqrt(3*dt);
        double nu   = r - div - 0.5*Math.pow(sigma,2);
        double edx  = Math.exp(dx);
        double xi   = Math.pow(sigma/dx,2);
        double pu   = -0.25*dt*(xi + nu/dx);
        double pm   = 1.0 + 0.5*dt*xi + 0.5*r*dt;
        double pd   = -0.25*dt*(xi - nu/dx);

        //Initialize asset prices at maturity
        double [] St = new double[2*Nj+1];
        St[0] = So*Math.exp(-N*dx);
        for (int j=1;j<=2*Nj;j++) {
            St[j] = St[j-1]*edx;
        }

        //Initialize option values at maturity
        double [][] V = new double[2][2*Nj+1];
        for (int j=0;j<=2*Nj;j++) {
            V[0][j] = poc.equals(Poc.CALL) ? Math.max(0.0,St[j]-K) : Math.max(0.0,K-St[j]);
        }

        //Compute derivative boundary conditions
        double lambdaL = poc.equals(Poc.CALL) ? 0.0 : -1 * (St[1]-St[0]);
        double lambdaU = poc.equals(Poc.CALL) ? (St[1]-St[0]) : 0.0;

        //Step backwards in time to now to get the initial option price
        for (int i=N-1; i>=0; i--) {
            solveCrankNicholsonTridiagonalSystem(V,pu,pm,pd,lambdaL,lambdaU);
            for(int j=0;j<=2*Nj;j++) {
                if (type.equals(OptionType.AMERICAN)) {
                    V[0][j] = poc.equals(Poc.CALL) ? Math.max(V[1][j],St[j]-K) : Math.max(V[1][j],K-St[j]);
                }
                else {
                    V[0][j] = V[1][j];
                }
            }
        }
        return new Valuation<Option>(option, V[0][Nj]);

    }

    private void solveCrankNicholsonTridiagonalSystem(double[][] V, double pu, double pm, double pd, double lambdaL, double lambdaU) {
        // substitute boundary condition at j=0 in j=1
        double [] pmp = new double[2*Nj];
        double [] pp  = new double[2*Nj];
        pmp[1] = pm + pd;
        pp[1]  = -pu*V[0][2] - (pm-2)*V[0][1] - pd*V[0][0] + pd*lambdaL;

        // eliminate upper diagonal
        for (int j=2;j<2*Nj;j++){
            pmp[j] = pm - pu*pd/pmp[j-1];
            pp[j]  = -pu*V[0][j+1] - (pm-2)*V[0][j] - pd*V[0][j-1] - pp[j-1]*pd/pmp[j-1];
        }

        // use boundary condition at j=2*Nj and equation at j=2*Nj-1
        V[1][2*Nj] = (pp[2*Nj-1] + pmp[2*Nj-1]*lambdaU )/(pu + pmp[2*Nj-1]);
        V[1][2*Nj-1] = V[1][2*Nj] - lambdaU;

        // back substitution
        for (int j=2*Nj-2;j>0;j--){
            V[1][j] = (pp[j]-pu*V[1][j+1])/pmp[j];
        }
        V[1][0] = V[1][1] - lambdaL; 

        return;
    }

    public Valuation getImpliedVolatility(Option instrument, double price) {
        return null;
    }
}