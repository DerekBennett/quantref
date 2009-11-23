import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.MathException;

import java.util.ArrayList;
import java.util.Random;
import java.util.TreeMap;

import umontreal.iro.lecuyer.probdist.ChiSquareNoncentralDist;
import com.sun.tools.javac.util.Pair;

/**
 * User: derekbennett
 * Date: Nov 8, 2009
 * Time: 1:29:33 PM
 */
public class hw3 {
    public static void main(String[] args) throws MathException {
        problem1a();
        problem1b();
        problem1c();

        problem2a();
        problem2b();
        problem2c();
        problem2d();

        problem3a();
        problem3b();
    }

    private static void problem1a() throws MathException {
        System.out.println("\n\nProblem 1.a");
        // Assume that the term structure of interest rate is ßat at 5% continuously compounded.
        // Assume that the forward bond price volatility is 8%.
        // Price a two year European call option with strike 0.75 on a seven-year pure discount bond.
        double testPrice = blackModelEurOptionOnPureDiscountBond(.05,.1,.8,1,5,0);
        System.out.printf("testPrice should be .0404 and actually is %4.4f\n",testPrice);
        double price = blackModelEurOptionOnPureDiscountBond(.05,.08,.75,2,7,0);
        System.out.printf("hw#3 Problem 1.a price = %4.4f\n",price);
    }

    private static double getPureDiscountBondPrice(double r, double s, double t) {
        return Math.exp(-r*(t-s));
    }
    private static double getForwardBondPrice(double r, double start, double optionExpiry, double bondMaturity) {
        double optionMaturityPrice =  getPureDiscountBondPrice(r,start,optionExpiry);
        double bondPrice = getPureDiscountBondPrice(r,start,bondMaturity);
        return bondPrice/optionMaturityPrice;
    }
    private static double blackModelEurOptionOnPureDiscountBond(double r,     // flat term structure interest rate
                                                                double sigma, // forward bond price volatility
                                                                double K,     // strike price
                                                                double T,     // time to option expiry
                                                                double s,     // time to bond maturity
                                                                double t)     // now
            throws MathException{
        double optionMaturityPrice = getPureDiscountBondPrice(r,t,T);
        double forwardBondPrice = getForwardBondPrice(r,t,T,s);
        double d1 = (Math.log(forwardBondPrice/K) + (Math.pow(sigma,2.0)/2.0)*(T-t))/(sigma*Math.sqrt(T-t));
        double d2 = d1 - sigma*Math.sqrt(T-t);
        NormalDistributionImpl N = new NormalDistributionImpl();
        return optionMaturityPrice*(forwardBondPrice*N.cumulativeProbability(d1) - K*N.cumulativeProbability(d2));
    }

    private static void problem1b() throws MathException {
        System.out.println("\n\nProblem 1.b");
        // Assume that the term structure of interest rate is ßat at 5% continuously compounded.
        // Price an interest rate cap at 4.5%. Assume that the cap is for
        // a two year period, the reset frequency is three months. Assume
        // that the volatility of the forward rate is 8%. The principal of the
        // swap is 1 million. (you may want to use software for this).

        // test caplet from p.191 Figure 6.5
        double testCapletPrice = blackModelCap(.05,.1,.045,.75,1,0,.25,1);
        System.out.printf("testPrice should be .0012 and actually is %4.4f\n\n",testCapletPrice);


        // each caplet is .25 in length due to reset frequency
        double r = .05;
        double TC = 2; //cap period
        double deltaTau = .25; // The length of each caplet
        double K = .045;// cap rate
        double sigma = .08; // volatility of forward rate
        double L = 1000000; // swap L
        double capPrice = 0.0;
        System.out.println("option   bond    forward    d1    d2    caplet");
        for (double T=deltaTau;T<TC;T+=deltaTau){
            double capletPrice = blackModelCap(r,sigma,K,T,T+deltaTau,0,deltaTau, L);
            capPrice+=capletPrice;
        }
        System.out.println("sum of caplets      = " + capPrice);
        System.out.println("Cap Price should be 8987.453 according to DerivaGem");
    }
    private static double blackModelCap(double r,     // flat term structure interest rate
                                        double sigma, // forward bond price volatility
                                        double K,     // strike price
                                        double T,     // time to option expiry
                                        double s,     // time to bond maturity
                                        double t,     // now
                                        double deltaTau,// reset period
                                        double L)     // Principal amount underlying the cap
            throws MathException{
        double optionMaturityPrice = getPureDiscountBondPrice(r,t,T);
        double bondMaturityPrice = getPureDiscountBondPrice(r,t,s);
        double forwardBondPrice = Math.log(optionMaturityPrice/bondMaturityPrice)/(s-T);//getForwardSwapRate(r,t,T,s);
        double d1 = (Math.log(forwardBondPrice/K) + (Math.pow(sigma,2.0)/2.0)*(T-t))/(sigma*Math.sqrt(T-t));
        double d2 = d1 - sigma*Math.sqrt(T-t);
        NormalDistributionImpl N = new NormalDistributionImpl();
        double caplet = bondMaturityPrice*(forwardBondPrice*N.cumulativeProbability(d1) - K*N.cumulativeProbability(d2))*deltaTau*L;
        System.out.printf("%1.4f    %1.4f    %1.4f    %1.4f    %1.4f    %4.3f\n",optionMaturityPrice,bondMaturityPrice,forwardBondPrice,d1,d2,caplet);
        return caplet;
    }
    private static double getForwardSwapRate(double r, double start, double optionExpiry, double bondMaturity) {
        double optionMaturityPrice =  getPureDiscountBondPrice(r,start,optionExpiry);
        double bondPrice = getPureDiscountBondPrice(r,start,bondMaturity);
        // terms in the log expression in inverse of the pure discount bond version of this method
        return Math.log(optionMaturityPrice/bondPrice)/(bondMaturity-optionExpiry);
    }


    private static void problem1c() throws MathException {
        System.out.println("\n\nProblem 1.c");
        // Assume that the term structure of interest rate is ßat at 5% continuously compounded.
        // Price a one year option which exercise into a new two year semi-annual payer swap.
        // The strike price of the option is 5% and the forward swap rate volatility is 15%.
        double testSwaption = blackModelSwaption(.05,.2,.05,2,1,0,.5);
        System.out.println("testSwaption = " + testSwaption);
        System.out.println("testSwaption should be .0052");

        double swaption = blackModelSwaption(.05,.15,.05,1,2,0,.5);
        System.out.println("swaption = " + swaption);
    }
    private static double blackModelSwaption(double r,     // flat term structure interest rate
                                             double sigma, // forward bond price volatility
                                             double K,     // strike price
                                             double T,     // time to option expiry
                                             double s,     // length of swap
                                             double t,     // now
                                             double deltaTau)// reset period of swap
            throws MathException
    {
        ArrayList<Double> allPDBs = new ArrayList<Double>();
        for (int i=0;i<=s/deltaTau;i++){
            double pdb = getPureDiscountBondPrice(r,0,T+i*deltaTau);
            allPDBs.add(pdb);
        }
        double totalForwardBondPrices = 0.0;
        double totalBondPrices = 0.0;
        for (int i=0;i<allPDBs.size()-1;i++){
            totalForwardBondPrices += allPDBs.get(i+1)/allPDBs.get(0);
            totalBondPrices+= allPDBs.get(i+1);
        }
        double Rfswap = (1-(allPDBs.get(allPDBs.size()-1)/allPDBs.get(0)));
        Rfswap = (Rfswap/totalForwardBondPrices)*2;
        double d1 = (Math.log(Rfswap/K) + (Math.pow(sigma,2.0)/2.0)*(T-t))/(sigma*Math.sqrt(T-t));
        double d2 = d1 - sigma*Math.sqrt(T-t);
        NormalDistributionImpl N = new NormalDistributionImpl();
        double swaption = deltaTau * totalBondPrices * (Rfswap*N.cumulativeProbability(d1) - K*N.cumulativeProbability(d2));
        return swaption;
    }


    private static void problem2a() throws MathException {
        System.out.println("\n\nProblem 2.a");
        // Price a two-year pure discount bond and a two-year European call option on a seven-year pure discount bond with strike price of 0.73.
        // Use the Vasicek model with r0 = 0.05, rbar = 0.05, alpha = 1.5. The volatility is sigma = 0.01.
        double testPrice = vasicekEurCallOptionOnPureDiscountBond(.15,.05,.01,.05,1,5,.67,0);
        System.out.printf("testPrice should be .1424 and actually is %4.4f\n",testPrice);

        double price = vasicekEurCallOptionOnPureDiscountBond(1.5,.05,.01,.05,2,7,.73,0);
        System.out.printf("problem price is %4.4f\n",price);
    }
    private static double vasicekEurCallOptionOnPureDiscountBond(double alpha,
                                             double rbar, // reversion mean
                                             double sigma,// volatility
                                             double r0,   // initial short rate
                                             double T,    // length of option on bond
                                             double s,    // length of PDB
                                             double K,    // strike
                                             double t)    // now    
            throws MathException{
        double alphaDiscount1 = 1-Math.exp(-alpha*(T-t));
        double sigma2 = Math.pow(sigma,2);
        double Bts = (1/alpha)*alphaDiscount1;
        double Rinf = rbar - .5*(sigma2/Math.pow(alpha,2));
        double lnAts = ((Rinf/alpha)*alphaDiscount1) - (T-t)*Rinf - (sigma2/(4*Math.pow(alpha,3)))*Math.pow(alphaDiscount1,2);
        double PtTPDB = Math.exp(lnAts)*Math.exp(-r0*Bts);
        System.out.println("PtTPDB = " + PtTPDB);

        double alphaDiscount2 = 1-Math.exp(-alpha*(s-t));
        sigma2 = Math.pow(sigma,2);
        Bts = (1/alpha)*alphaDiscount2;
        Rinf = rbar - .5*(sigma2/Math.pow(alpha,2));
        lnAts = ((Rinf/alpha)*alphaDiscount2) - (s-t)*Rinf - (sigma2/(4*Math.pow(alpha,3)))*Math.pow(alphaDiscount2,2);
        double PtsPDB = Math.exp(lnAts)*Math.exp(-r0*Bts);

        double sigmaR = (sigma/(alpha*(T-t)))*alphaDiscount1;

        double nutT = Math.sqrt((sigma2*(1-Math.exp(-2*alpha*(T-t))))/(2*alpha));
        double sigmaP = nutT*(1-Math.exp(-alpha*(s-T)))/alpha;

        double d1 = Math.log(PtsPDB/(K*PtTPDB))/sigmaP + sigmaP/2;
        double d2 = d1 - sigmaP;
        NormalDistributionImpl N = new NormalDistributionImpl();
        double price = PtsPDB*N.cumulativeProbability(d1) - K*PtTPDB*N.cumulativeProbability(d2);
        return price;
    }
    private static void problem2b() throws MathException {
        System.out.println("\n\nProblem 2.b");
        // Price a two-year pure discount bond and a two-year European call option on a seven-year pure discount bond with strike price of 0.73.
        // Use the CIR model with r0 = rbar = 0.05, alpha = 1.5, sigma = 0.1.
        double testPrice = cirEurCallOptionOnPureDiscountBond(.15,.05,.1,.05,1,5,.67,0);
        System.out.printf("testPrice should be .1463 and actually is %4.4f\n",testPrice);

        double price = cirEurCallOptionOnPureDiscountBond(1.5,.05,.1,.05,2,7,.73,0);
        System.out.printf("problem price is %4.4f\n",price);
    }
    private static double cirEurCallOptionOnPureDiscountBond(double alpha,
                                             double rbar, // reversion mean
                                             double sigma,// volatility
                                             double r0,   // initial short rate
                                             double T,    // length of option on bond
                                             double s,    // length of PDB
                                             double K,    // strike
                                             double t)    // now
            throws MathException{
        double sigma2 = Math.pow(sigma,2);
        double phi1=Math.sqrt(Math.pow(alpha,2) + (2*sigma2)); //.2062
        double phi2=(alpha+phi1)/2;                            //.1781
        double phi3=(2*alpha*rbar)/sigma2;                     //1.5

        double AtT = Math.pow( (phi1*Math.exp(phi2*(T-t))) / (phi2*(Math.exp(phi1*(T-t))-1) + phi1) ,phi3); //A(0,1) .9964
        double BtT = (Math.exp(phi1*(T-t))-1) / (phi2*(Math.exp(phi1*(T-t))-1) + phi1);                     //B(0,1) .9272
        double PtT = AtT * Math.exp(-BtT*r0);                                                               //P(0,1) .9513
        System.out.println("PtT = " + PtT);
        double Ats = Math.pow( (phi1*Math.exp(phi2*(s-t))) / (phi2*(Math.exp(phi1*(s-t))-1) + phi1) ,phi3);
        double Bts = (Math.exp(phi1*(s-t))-1) / (phi2*(Math.exp(phi1*(s-t))-1) + phi1);                     //B(0,5) 3.1499
        double Pts = Ats * Math.exp(-Bts*r0);                                                               //P(0,5) .7835

        double degreesOfFreedom=(4*alpha*rbar)/sigma2;                                                      // 3
        double theta=phi1;                                                                                  // theta .2602
        double phi4=(2*theta)/(sigma2*(Math.exp(theta*(T-t))-1));                                           // phi   180.09
        double psi=(alpha+theta)/sigma2;                                                                    // psi   35.62
        double ATs = Math.pow( (phi1*Math.exp(phi2*(s-T))) / (phi2*(Math.exp(phi1*(s-T))-1) + phi1) ,phi3); //A(1,5) .9521
        double BTs = (Math.exp(phi1*(s-T))-1) / (phi2*(Math.exp(phi1*(s-T))-1) + phi1);                     //B(1,5) 2.9498
        double nonCentrality1=(2*Math.pow(phi4,2)*r0*Math.exp(theta*(T-t)))/(phi4+psi+BTs);                 //       18.2288
        double nonCentrality2=(2*Math.pow(phi4,2)*r0*Math.exp(theta*(T-t)))/(phi4+psi);                     //       18.4781
        double rstar=Math.log(ATs/K)/BTs;                                                                   // r*    .1191
        double d1=2*rstar*(phi4+psi+BTs);                                                                   //       52.0917
        double d2=2*rstar*(phi4+psi);                                                                       //       51.3890
        double Chi21 = ChiSquareNoncentralDist.cdf(degreesOfFreedom,nonCentrality1,d1);                     //       .9999395
        double Chi22 = ChiSquareNoncentralDist.cdf(degreesOfFreedom,nonCentrality2,d2);                     //       .9999059
        double price = Pts*Chi21 - K*PtT*Chi22;
        return price;
    }
    private static void problem2c() throws MathException {
        System.out.println("\n\nProblem 2.c");
        // Price a two-year pure discount bond and a two-year European call option on a seven-year pure discount bond with strike price of 0.73.
        // Using the Ho-Lee model with r(1) = 0.05 and ? = 0.01 calculate the discount bond price above after 1 year, and the call option described above.
        double testPrice = hoLeeEurCallOptionOnPureDiscountBond(.05,.01,.05,1,5,.8187,0);
        double testBondPrice = hoLeePureDiscountBond(.05,.01,.05,1,5,0);
        System.out.printf("testPrice should be     .0124 and actually is %4.4f\n",testPrice);
        System.out.printf("testBondPrice should be .8181 and actually is %4.4f\n",testBondPrice);

        double price = hoLeeEurCallOptionOnPureDiscountBond(.05,.01,.05,2,7,.73,0);
        double bondPrice = hoLeePureDiscountBond(.05,.01,.05,2,7,0);
        System.out.printf("problem price is           %4.4f\n",price);
        System.out.printf("problem bond price(t=0) is %4.4f\n",bondPrice);
        bondPrice = hoLeePureDiscountBond(.05,.01,.05,2,7,1);
        System.out.printf("problem bond price(t=1) is %4.4f\n",bondPrice);
    }

    private static double hoLeePureDiscountBond(double r, // flat term-structure rate
                                             double sigma,// volatility
                                             double r1,   // short rate after one year
                                             double T,    // length of option on bond
                                             double s,    // length of PDB
                                             double t)    // now
            throws MathException{
        double PtT = Math.exp(-r*T);
        double Pts = Math.exp(-r*s);
        double deltaT = .1;
        double BTs = (s-T);
        double PtTPlusDeltaT = Math.exp(-r*(T+deltaT));
        double PtTMinusDeltaT = Math.exp(-r*(T-deltaT));
        double slope = (Math.log(PtTPlusDeltaT) - Math.log(PtTMinusDeltaT))/(2*deltaT);
        double lnATs = Math.log(Pts/PtT) - BTs*(slope) - .5*Math.pow(sigma,2)*(T-t)*Math.pow(BTs,2);
        double price = Math.exp(lnATs)*Math.exp(-BTs*r1);
        return price;
    }
    private static double hoLeeEurCallOptionOnPureDiscountBond(double r, // flat term-structure rate
                                             double sigma,// volatility
                                             double r1,   // short rate after one year
                                             double T,    // length of option on bond
                                             double s,    // length of PDB
                                             double K,    // strike
                                             double t)    // now
            throws MathException{
        double PtT = Math.exp(-r*T);
        double Pts = Math.exp(-r*s);
        double sigmaP = sigma*(s-T)*Math.sqrt(T-t);
        double d1 =  Math.log(Pts/(K*PtT))/(sigmaP) + (sigmaP/2);
        double d2 =  d1 - sigmaP;
        NormalDistributionImpl N = new NormalDistributionImpl();
        double price = Pts*N.cumulativeProbability(d1) - K*PtT*N.cumulativeProbability(d2);
        return price;
    }

    private static void problem2d() throws MathException {
        System.out.println("\n\nProblem 2.d");
        // Price a two-year pure discount bond and a two-year European call option on a seven-year pure discount bond with strike price of 0.73.
        // Use the Hull-White model with rbar = 0.05, ? = .1, and ? = 0.01.
        // Calculate the discount bond price after 1 year, and the call option as in part c).
        double testPrice = hullWhiteEurCallOptionOnPureDiscountBond(.05,.1,.01,.05,1,5,.8187,0);
        System.out.printf("testPrice should be     .0098 and actually is %4.4f\n",testPrice);
        double testBondPrice = hullWhitePureDiscountBond(.05,.1,.01,.05,1,5,0);
        System.out.printf("testBondPrice should be .8183 and actually is %4.4f\n",testBondPrice);

        double price = hullWhiteEurCallOptionOnPureDiscountBond(.05,.1,.01,.05,2,7,.73,0);
        System.out.printf("problem price is           %4.4f\n",price);
        double bondPrice = hullWhitePureDiscountBond(.05,.1,.01,.05,2,7,0);
        System.out.printf("problem bond price(t=0) is %4.4f\n",bondPrice);
        bondPrice = hullWhitePureDiscountBond(.05,.1,.01,.05,2,7,1);
        System.out.printf("problem bond price(t=1) is %4.4f\n",bondPrice);
    }

    private static double hullWhitePureDiscountBond(double r, // flat term-structure rate
                                             double alpha,// inverse mean-reversion speed parameter
                                             double sigma,// volatility
                                             double r1,   // short rate after one year
                                             double T,    // length of option on bond
                                             double s,    // length of PDB
                                             double t)    // now
            throws MathException{
        double PtT = Math.exp(-r*T);
        double Pts = Math.exp(-r*s);
        double deltaT = .1;
        double BTs = (1/alpha)*(1-Math.exp(-alpha*(s-T)));
        double PtTPlusDeltaT = Math.exp(-r*(T+deltaT));
        double PtTMinusDeltaT = Math.exp(-r*(T-deltaT));
        double slope = (Math.log(PtTPlusDeltaT) - Math.log(PtTMinusDeltaT))/(2*deltaT);
        double lnATs = Math.log(Pts/PtT) - BTs*(slope) - (1/(4*Math.pow(alpha,3)))*Math.pow(sigma,2)*Math.pow((Math.exp(-alpha*(s-t))-Math.exp(-alpha*(T-t))),2)*(Math.exp(2*alpha*(T-t))-1);
        double price = Math.exp(lnATs)*Math.exp(-BTs*r1);
        return price;
    }
    private static double hullWhiteEurCallOptionOnPureDiscountBond(double r, // flat term-structure rate
                                             double alpha,// inverse mean-reversion speed parameter 
                                             double sigma,// volatility
                                             double r1,   // short rate after one year
                                             double T,    // length of option on bond
                                             double s,    // length of PDB
                                             double K,    // strike
                                             double t)    // now
            throws MathException{
        double PtT = Math.exp(-r*T);
        System.out.println("PtT = " + PtT);
        double Pts = Math.exp(-r*s);
        double sigmaP = Math.sqrt((Math.pow(sigma,2)/(2*Math.pow(alpha,3))) * (1-Math.exp(-2*alpha*(T-t))) * Math.pow(1-Math.exp(-alpha*(s-T)),2));
        double d1 =  Math.log(Pts/(K*PtT))/(sigmaP) + (sigmaP/2);
        double d2 =  d1 - sigmaP;
        NormalDistributionImpl N = new NormalDistributionImpl();
        double price = Pts*N.cumulativeProbability(d1) - K*PtT*N.cumulativeProbability(d2);
        return price;
    }

    private static void problem3a() throws MathException {
        System.out.println("\n\nProblem 3.a");
        /*
            Price a one year option which exercise into a three year semi-annual payer swap
            (note that unlike in part 1b) the swap is going on while you hold the option).
            The strike price of the option is 5% and r0 = 0.05.
            Simulate the interest rate paths using the CIR model,
            and calculate the swaption value as the discounted expectation of future payments.

            CIR model params: r0 = rbar = 0.05, ? = 1.5, ? = 0.1.
         */
        long start=System.currentTimeMillis();
        double T=1;
        int s=3;
        double couponPeriod=.5;
        double r0=.05,rbar=.05,K=.05;
        double alpha=1.5;
        double sigma=.1;

        // Build N paths of M samples each
        int N=0;
        int M=1000*s;//If M is a multiple of the number of coupon dates, then it will be simple to pick them out
        double dt=s*1.0/M;
        // each path contains a map of (coupon date,rate) tuples
        ArrayList<TreeMap<Double,Double>> paths= new ArrayList<TreeMap<Double, Double>>();
        Random random = new Random();
        for (int i=0; i<=N;i++){
            // Record interest rates at each of the coupon dates
            double r=r0;
            TreeMap<Double,Double> dateRates = new TreeMap<Double, Double>();
            paths.add(dateRates);
            for(int j=0;j<=M;j++){
                // Use Euler's discretization of the CIR "square-root" process
                double rPrime = r + (alpha * (rbar - r) * dt) + (sigma * Math.sqrt(r) * Math.sqrt(dt) * random.nextGaussian());
                r = Math.max(0,rPrime);// Prevent unreal rates
                System.out.printf("%f\t%f\n",j*dt,r);
                if (j*dt % .5 == 0){
                    //System.out.printf("t=%f\tr = %f\n",j*dt,r);
                    dateRates.put(j*dt,r);
                }
            }
        }
        //Now, for each path, compute the option and swap values with a running average
        double accumulatedSwaptionValue=0.0;
        for(TreeMap<Double, Double> path : paths){
            ArrayList<Double> allPDBs = new ArrayList<Double>();
            //Step through each coupon
            double swaptionValue=0.0;
            for(Double couponDate:path.keySet()){
                // Calculate the expected payment at that date. Discount it back to the present and add it to the total
                double pdb = getPureDiscountBondPrice(path.get(couponDate),0,couponDate);
                allPDBs.add(pdb);
            }
            double totalForwardBondPrices = 0.0;
            double totalBondPrices = 0.0;
            for (int i=0;i<allPDBs.size()-1;i++){
                totalForwardBondPrices += allPDBs.get(i+1)/allPDBs.get(0);
                totalBondPrices+= allPDBs.get(i+1);
            }
            double Rfswap = (1-(allPDBs.get(allPDBs.size()-1)/allPDBs.get(0)));
            Rfswap = (Rfswap/totalForwardBondPrices)*2;
            double d1 = (Math.log(Rfswap/K) + (Math.pow(sigma,2.0)/2.0)*(T-0))/(sigma*Math.sqrt(T-0));
            double d2 = d1 - sigma*Math.sqrt(T-0);
            NormalDistributionImpl NormDist = new NormalDistributionImpl();
            swaptionValue = couponPeriod * totalBondPrices * (Rfswap*NormDist.cumulativeProbability(d1) - K*NormDist.cumulativeProbability(d2));
//            System.out.println("swaptionValue = " + swaptionValue);
            accumulatedSwaptionValue+=swaptionValue;
        }
        System.out.println("M = " + M);
        System.out.println("N = " + N);
        System.out.println("average swaptionValue = " + accumulatedSwaptionValue/paths.size());
        System.out.println("calculation took " + (System.currentTimeMillis()-start) + "ms");
    }
    private static void problem3b() throws MathException {
        System.out.println("\n\nProblem 3.b");
        /*
            Price a one year option which exercise into a three year semi-annual payer swap
            (note that unlike in part 1b) the swap is going on while you hold the option).
            The strike price of the option is 5% and r0 = 0.05.
            Simulate the interest rate paths using the Longsta? and Schwartz model,
            and calculate the swaption value as the discounted expectation of future payments.

            Model params: ? = 0.0024, ? = 0.0655,  ? = 0.24581, ? = 0.0164, ? = 0.24537, ? = 1.1457, r0 = rbar = 0.05, ? = 0.001.
         */
        long start=System.currentTimeMillis();
        double T=1;
        int s=3;
        double couponPeriod=.5;
        double K    = .05;
        double alpha= .0024;
        double beta = .0655;
        double gamma= .24581;
        double delta= 0.0164;
        double eta  = 0.24537;
        double theta= 1.1457;
        double r0   = 0.05;
        double rbar = 0.05;
        double nu   = 0.001;

        // Build N paths of M samples each
        int N=0;
        int M=1000*s;//If M is a multiple of the number of coupon dates, then it will be simple to pick them out
        double dt=s*1.0/M;
        // each path contains a map of (coupon date,rate) tuples
        ArrayList<TreeMap<Double,Double>> paths= new ArrayList<TreeMap<Double,Double>>();
        Random random = new Random();
        for (int i=0; i<=N;i++){
            // Record interest rates at each of the coupon dates
            double x=1;
            double y=(r0-alpha*x)/beta;
            TreeMap<Double,Double> dateRates = new TreeMap<Double,Double>();
            paths.add(dateRates);
            for(int j=0;j<=M;j++){
                // Use Euler's discretization of the Longstaff and Schwartz two-factor stochastic volatility process
                double xprime = x + (gamma-delta*x)*dt + Math.sqrt(x)*Math.sqrt(dt)*random.nextGaussian();
                double yprime = y + (eta-theta*y)*dt + Math.sqrt(y)*Math.sqrt(dt)*random.nextGaussian();
                x=xprime;
                y=yprime;
                double r = alpha*x + beta*y;
                if (j*dt % .5 == 0){
                    System.out.printf("t=%f\tr = %f\n",j*dt,r);
                    dateRates.put(j*dt,r);
                }
            }
        }
        //Now, for each path, compute the option and swap values with a running average
        double accumulatedSwaptionValue=0.0;
        for(TreeMap<Double, Double> path : paths){
            ArrayList<Double> allPDBs = new ArrayList<Double>();
            //Step through each coupon
            double swaptionValue=0.0;
            for(Double couponDate:path.keySet()){
                // Calculate the expected payment at that date. Discount it back to the present and add it to the total
                double pdb = getPureDiscountBondPrice(path.get(couponDate),0,couponDate);
                allPDBs.add(pdb);
            }
            double totalForwardBondPrices = 0.0;
            double totalBondPrices = 0.0;
            for (int i=0;i<allPDBs.size()-1;i++){
                totalForwardBondPrices += allPDBs.get(i+1)/allPDBs.get(0);
                totalBondPrices+= allPDBs.get(i+1);
            }
            double Rfswap = (1-(allPDBs.get(allPDBs.size()-1)/allPDBs.get(0)));
            Rfswap = (Rfswap/totalForwardBondPrices)*2;
            double d1 = (Math.log(Rfswap/K) + (Math.pow(nu,2.0)/2.0)*(T-0))/(nu*Math.sqrt(T-0));
            double d2 = d1 - nu*Math.sqrt(T-0);
            NormalDistributionImpl NormDist = new NormalDistributionImpl();
            swaptionValue = couponPeriod * totalBondPrices * (Rfswap*NormDist.cumulativeProbability(d1) - K*NormDist.cumulativeProbability(d2));
//            System.out.println("swaptionValue = " + swaptionValue);
            accumulatedSwaptionValue+=swaptionValue;
        }
        System.out.println("M = " + M);
        System.out.println("N = " + N);
        System.out.println("average swaptionValue = " + accumulatedSwaptionValue/paths.size());
        System.out.println("calculation took " + (System.currentTimeMillis()-start) + "ms");
    }
}
