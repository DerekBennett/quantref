import quantref.*;
import static quantref.OptionType.*;
import static quantref.Poc.*;

import java.util.*;

import org.joda.time.LocalDate;

/**
 * User: derekbennett
 * Date: Oct 18, 2009
 * Time: 2:36:26 PM
 */
public class hw2 {
    public static void main(String[] args) throws InvalidInstrumentException, ImmutablesAreImmutableException {
        unitTests();
        hw2_problem1();
        hw2_problem2();
        hw2_problem3();
        hw2_problem4();
        hw2_problem6();
        hw2_problem7();
        hw2_problem9();
        hw2_problem12and13();
        hw2_problem14();
    }

    private static void hw2_problem1() {
        System.out.println("\n\n===================================================");
        System.out.println("PROBLEM 1");
        System.out.println("Read Ch. 3 in Clewlow. Implement a trinomial tree scheme to calculate \n" +
                "values for call and put options both European and American.\n");
        System.out.println("ANSWER");
        System.out.println("Please refer to TrinomialTreeValuator.java in the submitted code files\n");
    }

    private static void hw2_problem2() throws InvalidInstrumentException, ImmutablesAreImmutableException {
        Linear syntheticAmPutUnderlying = new Linear("syntheticAmPutUnderlying",null,0.0);
        LocalDate TODAY = new LocalDate(2009,1,1);
        LocalDate EXPIRY = TODAY.plusMonths(1); // This makes T=1/12
        Option synthPutAmOption = new Option(PUT, AMERICAN, TODAY,102,EXPIRY,syntheticAmPutUnderlying); //K=102
        Option synthPutEurOption = new Option(PUT, EUROPEAN, TODAY,102,EXPIRY,syntheticAmPutUnderlying); //K=102
        BaseStateOfTheWorld sotw = new BaseStateOfTheWorld(TODAY,.01); //r=.01
        HashMap<Instrument, Double> prices = new HashMap<Instrument, Double>();
        HashMap<Option, Double> vols = new HashMap<Option, Double>();
        prices.put(syntheticAmPutUnderlying,100.0); //S=100
        vols.put(synthPutAmOption,.25); //sigma=.25
        vols.put(synthPutEurOption,.25); //sigma=.25
        sotw.setPrices(prices);
        sotw.setVolatilities(vols);

        BlackScholesValuator bsv = new BlackScholesValuator();
        bsv.setStateOfTheWorld(sotw);
        Valuation bsVal = bsv.getValuation(synthPutEurOption);

        BinomialTreeValuator btv = new BinomialTreeValuator(sotw,200);
        Valuation btVal200 = btv.getValuation(synthPutAmOption);

        TrinomialTreeValuator ttv = new TrinomialTreeValuator(sotw,200);
        Valuation ttVal200 = ttv.getValuation(synthPutAmOption);

        System.out.println("\n\n===================================================");
        System.out.println("PROBLEM 2");
        System.out.println("Refer to the Þrst assignment. Use a synthetic American Put option with \n" +
                "r = 0.01, T = 1/12, ? = 0.25, S0 = 100, K = 102. Treat this option as \n" +
                "European and calculate its price using the Black-Scholes formula. Then \n" +
                "treat it properly as an American option and price this option using the \n" +
                "binomial tree you implemented in the last assignment. Then price the \n" +
                "option using the Trinomial tree from the previous part. Compare the \n" +
                "three prices you obtained.\n");
        System.out.println("ANSWER");
        System.out.println(" Values of a synthetic American Put option with \n" +
                           "r = 0.01, T = 1/12, sigma = 0.25, S0 = 100, K = 102\n");
        System.out.println("Black-Scholes  price       = " + bsVal.getPv());
        System.out.println("Binomial Tree  price@200   = " + btVal200.getPv());
        System.out.println("Trinomial Tree price@200   = " + ttVal200.getPv());
        System.out.println("\n\n");
    }

    private static void hw2_problem3() throws ImmutablesAreImmutableException, InvalidInstrumentException {
        Linear syntheticAmPutUnderlying = new Linear("syntheticAmPutUnderlying",null,0.0);
        LocalDate TODAY = new LocalDate(2009,1,1);
        LocalDate EXPIRY = TODAY.plusMonths(1); // This makes T=1/12
        Option synthPutAmOption = new Option(PUT, AMERICAN, TODAY,102,EXPIRY,syntheticAmPutUnderlying); //K=102
        BaseStateOfTheWorld sotw = new BaseStateOfTheWorld(TODAY,.01); //r=.01
        HashMap<Instrument, Double> prices = new HashMap<Instrument, Double>();
        HashMap<Option, Double> vols = new HashMap<Option, Double>();
        prices.put(syntheticAmPutUnderlying,100.0); //S=100
        vols.put(synthPutAmOption,.25); //sigma=.25
        sotw.setPrices(prices);
        sotw.setVolatilities(vols);

        System.out.println("\n\n===================================================");
        System.out.println("PROBLEM 3");
        System.out.println("Again refer to part 2. Now calculate the American put option price using \n" +
                "the binomial tree with 50, 51, 52 and 53 time steps. What do you see? \n" +
                "Then Þnd out how many steps are needed for the trinomial tree to Þnd \n" +
                "answers of the same order of magnitude as the binomial tree with 50 steps. \n" +
                "(Try di?erent numbers of time steps until you Þnd a put value similar with \n" +
                "the one obtained for the binomial tree).\n");
        System.out.println("ANSWER\n");
        for (int i=50;i<54;i++){
            BinomialTreeValuator btv = new BinomialTreeValuator(sotw,i);
            Valuation btVal = btv.getValuation(synthPutAmOption);
            System.out.println("binomial value with " + i +" steps = " + btVal.getPv());
        }
        System.out.println("The binomial values seem to be oscillating about the 3.99... convergence.");
        System.out.println("\n");
        System.out.println("Trinomial values with ");
        for (int i=4;i<101;i++){
            TrinomialTreeValuator ttv = new TrinomialTreeValuator(sotw,i);
            Valuation ttVal = ttv.getValuation(synthPutAmOption);
            System.out.println(i +" steps = " + ttVal.getPv() +
                    (Math.abs(ttVal.getPv()-3.9865580663557982)<.001 ? "\t value close to value of Binomial at 50 steps" : ""));
        }
        System.out.println("\nBecause of the greater number of branches (and therefore pricing paths), \n" +
                "there will be a quicker convergence using the Trinomial tree. The Trinomial reaches \n" +
                "a similar value to that of the Binomial in 32 steps versus the Binomial's 50 steps. \n" +
                "Further, as the behavior of both models is cyclic and convergent, either model might \n" +
                "have a number closer to the convergence value at any given number of steps.");
    }
    
    private static void hw2_problem4() throws InvalidInstrumentException, ImmutablesAreImmutableException {
        Linear testEurUnderlying = new Linear("syntheticUnderlying",null,0.03);
        LocalDate TODAY = new LocalDate(2009,1,1);
        LocalDate EXPIRY = TODAY.plusMonths(12); //T=1
        Option testEurCallOption = new Option(CALL, EUROPEAN, TODAY,100,EXPIRY,testEurUnderlying);
        Option testAmPutOption = new Option(PUT, AMERICAN, TODAY,100,EXPIRY,testEurUnderlying);
        BaseStateOfTheWorld sotw = new BaseStateOfTheWorld(TODAY,.06);
        HashMap<Instrument, Double> prices = new HashMap<Instrument, Double>();
        HashMap<Option, Double> vols = new HashMap<Option, Double>();
        prices.put(testEurUnderlying,100.0);
        vols.put(testEurCallOption,.2);
        prices.put(testAmPutOption,100.0);
        vols.put(testAmPutOption,.2);
        sotw.setPrices(prices);
        sotw.setVolatilities(vols);

        ExplicitFiniteDifferenceValuator efdv = new ExplicitFiniteDifferenceValuator(sotw,3,3);
        Valuation efdValCall = efdv.getValuation(testEurCallOption);
        Valuation efdValPut = efdv.getValuation(testAmPutOption);
        System.out.println("\n\n===================================================");
        System.out.println("PROBLEM 4");
        System.out.println("Implement an explicit Þnite difference method for the same purpose as the \n" +
                "part (1).\n");
        System.out.println("ANSWER\n");
        System.out.println("Please see ExplicitFiniteDifferenceValuator.java and the hw2_problem4 in hw2.java");
        System.out.println("European call example value from p.61 is " + efdValCall.getPv() + " and should be 8.5455");
        System.out.println("American put example  value from p.63 is " + efdValPut.getPv() + " and should be 6.0058");
    }

    private static void hw2_problem6() throws InvalidInstrumentException, ImmutablesAreImmutableException {
        //SP500 CALL OCT09@1070 (market price=(22.7,25.1),implied vol=0.195313)
        System.out.println("\n\n===================================================");
        System.out.println("PROBLEM 6");
        Linear SP500 = new Linear("SP500",null,0.0);
        LocalDate OCT09_DMO = new LocalDate(2009,10,1);
        LocalDate OCT09_EXP_SP500 = new LocalDate(2009,10,17);
        Option SP500Option = new Option(CALL,EUROPEAN,OCT09_DMO,1070,OCT09_EXP_SP500,SP500);//ATM option
        BaseStateOfTheWorld sotw = new BaseStateOfTheWorld(new LocalDate(2009,9,16),0.0017);// Fed Funds rate of 0.17%
        HashMap<Instrument, Double> prices = new HashMap<Instrument, Double>();
        prices.put(SP500,1068.76);
        sotw.setPrices(prices);
        HashMap<Option, Double> vols = new HashMap<Option, Double>();
        vols.put(SP500Option,0.1953125);
        sotw.setVolatilities(vols);

        BlackScholesValuator bsv = new BlackScholesValuator();
        bsv.setStateOfTheWorld(sotw);
        double bsSP500Price=bsv.getValuation(SP500Option).getPv();
        System.out.println("bsSP500Price = " + bsSP500Price);
        BinomialTreeValuator btv;
        TrinomialTreeValuator ttv;
        ExplicitFiniteDifferenceValuator efdv;
        System.out.println("steps\tbinomial\ttrinomial\texplicitFD");
        for (int steps=10; steps<401; steps++){
            btv = new BinomialTreeValuator(sotw,steps);
            ttv = new TrinomialTreeValuator(sotw,steps);
            efdv = new ExplicitFiniteDifferenceValuator(sotw,steps,steps);
            double btvSP500Diff  = btv.getValuation(SP500Option).getPv();
            double ttvSP500Diff  = ttv.getValuation(SP500Option).getPv();
            double efdvSP500Diff = efdv.getValuation(SP500Option).getPv();
            System.out.printf("%d\t%f\t%f\t%f\n",steps,btvSP500Diff,ttvSP500Diff,efdvSP500Diff);
        }
    }

    private static void hw2_problem7() throws InvalidInstrumentException, ImmutablesAreImmutableException {
        Linear testEurUnderlying = new Linear("syntheticUnderlying",null,0.0);
        LocalDate TODAY = new LocalDate(2009,1,1);
        LocalDate EXPIRY = TODAY.plusMonths(12); //T=1
        Option testEurCallOption = new Option(CALL, EUROPEAN, TODAY,100,EXPIRY,testEurUnderlying);
        BaseStateOfTheWorld sotw = new BaseStateOfTheWorld(TODAY,.06);
        HashMap<Instrument, Double> prices = new HashMap<Instrument, Double>();
        HashMap<Option, Double> vols = new HashMap<Option, Double>();
        prices.put(testEurUnderlying,100.0);
        vols.put(testEurCallOption,.2);
        sotw.setPrices(prices);
        sotw.setVolatilities(vols);

        System.out.println("\n\n===================================================");
        System.out.println("PROBLEM 7");
        System.out.println("Implement a simple Monte Carlo scheme using n simulation trials for \n" +
                "European Put and Call options. Implement a calculation of the standard \n" +
                "error of the estimate (use 300 time intervals).\n");
        System.out.println("ANSWER\n");
        System.out.println("Please see MonteCarloValuator.java and the hw2_problem7 in hw2.java");

        //Make the i limit equal to 10000000 to see the convergence
        for (int i=100; i<=10000;i=i*10){
            MonteCarloValuator mcv = new MonteCarloValuator(sotw,300,i);
            long start = System.currentTimeMillis();
            Valuation mcvValCall = mcv.getValuation(testEurCallOption);
            long total =  System.currentTimeMillis()-start;
            System.out.println("number of iterations = " + i + " took " + total + "ms");
            System.out.println("Euro call is " + mcvValCall.getPv() + " and should be 10.98954699");
            System.out.println("mcv.stddev = " + mcv.stddev);
            System.out.println("mcv.stderr = " + mcv.stderr);
            System.out.println("Estimate of number of trials to get stderr=.001: " +Math.pow(mcv.stddev/.001,2));
        }
    }

    private static void hw2_problem9() throws InvalidInstrumentException, ImmutablesAreImmutableException {
        System.out.println("\n\n===================================================");
        System.out.println("PROBLEM 9");
        // Underlying stock objects
        Linear GOOG = new Linear("GOOG",null,0.0);
        Linear SP500 = new Linear("SP500",null,0.0);

        // Important dates
        LocalDate OCT09_DMO = new LocalDate(2009,10,1);
        LocalDate DEC09_DMO = new LocalDate(2009,12,1);
        LocalDate OCT09_EXP_GOOG = new LocalDate(2009,10,16);
        LocalDate OCT09_EXP_SP500 = new LocalDate(2009,10,17);
        LocalDate DEC09_EXP_GOOG = new LocalDate(2009,12,18);
        LocalDate DEC09_EXP_SP500 = new LocalDate(2009,12,19);

        //GOOG ATM options
        Option GOOG_OCT_CALL_Option = new Option(CALL,EUROPEAN,OCT09_DMO,490,OCT09_EXP_GOOG,GOOG);
        Option GOOG_OCT_PUT_Option  = new Option(PUT,EUROPEAN,OCT09_DMO,490,OCT09_EXP_GOOG,GOOG);
        Option GOOG_DEC_CALL_Option = new Option(CALL,EUROPEAN,DEC09_DMO,490,DEC09_EXP_GOOG,GOOG);
        Option GOOG_DEC_PUT_Option  = new Option(PUT,EUROPEAN,DEC09_DMO,490,DEC09_EXP_GOOG,GOOG);
        //SP500 ATM options
        Option SP500_OCT_CALL_Option = new Option(CALL,EUROPEAN,OCT09_DMO,1070,OCT09_EXP_SP500,SP500);
        Option SP500_OCT_PUT_Option  = new Option(PUT,EUROPEAN,OCT09_DMO,1070,OCT09_EXP_SP500,SP500);
        Option SP500_DEC_CALL_Option = new Option(CALL,EUROPEAN,DEC09_DMO,1070,DEC09_EXP_SP500,SP500);
        Option SP500_DEC_PUT_Option  = new Option(PUT,EUROPEAN,DEC09_DMO,1070,DEC09_EXP_SP500,SP500);

        //Set-up the state-of-the-world (along with implied volatilities because my calculations from hw1 were incorrect)
        BaseStateOfTheWorld sotw = new BaseStateOfTheWorld(new LocalDate(2009,9,16),0.0017);// Fed Funds rate of 0.17%
        BlackScholesValuator bsv = new BlackScholesValuator();
        bsv.setStateOfTheWorld(sotw);
        HashMap<Instrument, Double> prices = new HashMap<Instrument, Double>();
        prices.put(GOOG,488.29);
        prices.put(SP500,1068.76);
        sotw.setPrices(prices);
        HashMap<Option, Double> vols = new HashMap<Option, Double>();
        vols.put(GOOG_OCT_CALL_Option ,bsv.getImpliedVolatility(GOOG_OCT_CALL_Option ,(15.8+16.1)/2.0).getPv());
        vols.put(GOOG_OCT_PUT_Option  ,bsv.getImpliedVolatility(GOOG_OCT_PUT_Option  ,(17.3+17.6)/2.0).getPv());
        vols.put(GOOG_DEC_CALL_Option ,bsv.getImpliedVolatility(GOOG_DEC_CALL_Option ,(26.1+26.6)/2.0).getPv());
        vols.put(GOOG_DEC_PUT_Option  ,bsv.getImpliedVolatility(GOOG_DEC_PUT_Option  ,(27.5+27.9)/2.0).getPv());

        vols.put(SP500_OCT_CALL_Option,bsv.getImpliedVolatility(SP500_OCT_CALL_Option,(22.7+25.1)/2.0).getPv());
        vols.put(SP500_OCT_PUT_Option ,bsv.getImpliedVolatility(SP500_OCT_PUT_Option ,(25.8+28.2)/2.0).getPv());
        vols.put(SP500_DEC_CALL_Option,bsv.getImpliedVolatility(SP500_DEC_CALL_Option,(43.7+46.1)/2.0).getPv());
        vols.put(SP500_DEC_PUT_Option ,bsv.getImpliedVolatility(SP500_DEC_PUT_Option ,(50.1+52.5)/2.0).getPv());
        sotw.setVolatilities(vols);

        //Set-up a portfolio to iterate over
        List<Option> portfolio = new ArrayList<Option>();
        portfolio.add(GOOG_OCT_CALL_Option );
        portfolio.add(GOOG_OCT_PUT_Option  );
        portfolio.add(GOOG_DEC_CALL_Option );
        portfolio.add(GOOG_DEC_PUT_Option  );
        portfolio.add(SP500_OCT_CALL_Option);
        portfolio.add(SP500_OCT_PUT_Option );
        portfolio.add(SP500_DEC_CALL_Option);
        portfolio.add(SP500_DEC_PUT_Option );

        int timeSteps=300;
        BinomialTreeValuator btv = new BinomialTreeValuator(sotw,timeSteps);
        TrinomialTreeValuator ttv= new TrinomialTreeValuator(sotw,timeSteps);
        ExplicitFiniteDifferenceValuator efdv= new ExplicitFiniteDifferenceValuator(sotw,timeSteps,timeSteps);
        MonteCarloValuator mcv= new MonteCarloValuator(sotw,timeSteps,100000);
        System.out.println("OPTION\tBlack-Scholes\tBinomial\tTrinomial\tEFD\tMonteCarlo");
        for (Option option : portfolio){
            System.out.printf("%s\t%f\t%f\t%f\t%f\t%f\n",option.getDescription(),
                    bsv.getValuation(option).getPv(),
                    btv.getValuation(option).getPv(),
                    ttv.getValuation(option).getPv(),
                    efdv.getValuation(option).getPv(),
                    mcv.getValuation(option).getPv());
        }
    }

    private static void hw2_problem12and13() throws InvalidInstrumentException, ImmutablesAreImmutableException {
        System.out.println("\n\n===================================================");
        System.out.println("PROBLEM 12 and 13");
        // Underlying stock objects
        Linear GOOG = new Linear("GOOG",null,0.0);
        Linear SP500 = new Linear("SP500",null,0.0);

        // Important dates
        LocalDate OCT09_DMO = new LocalDate(2009,10,1);
        LocalDate DEC09_DMO = new LocalDate(2009,12,1);
        LocalDate OCT09_EXP_GOOG = new LocalDate(2009,10,16);
        LocalDate OCT09_EXP_SP500 = new LocalDate(2009,10,17);
        LocalDate DEC09_EXP_GOOG = new LocalDate(2009,12,18);
        LocalDate DEC09_EXP_SP500 = new LocalDate(2009,12,19);

        //GOOG ATM options
        Option GOOG_OCT_PUT_Option  = new Option(PUT,AMERICAN,OCT09_DMO,490,OCT09_EXP_GOOG,GOOG);
        Option GOOG_DEC_PUT_Option  = new Option(PUT,AMERICAN,DEC09_DMO,490,DEC09_EXP_GOOG,GOOG);
        //SP500 ATM options
        Option SP500_OCT_PUT_Option  = new Option(PUT,AMERICAN,OCT09_DMO,1070,OCT09_EXP_SP500,SP500);
        Option SP500_DEC_PUT_Option  = new Option(PUT,AMERICAN,DEC09_DMO,1070,DEC09_EXP_SP500,SP500);

        //Set-up the state-of-the-world (along with implied volatilities because my calculations from hw1 were incorrect)
        BaseStateOfTheWorld sotw = new BaseStateOfTheWorld(new LocalDate(2009,9,16),0.0017);// Fed Funds rate of 0.17%
        BlackScholesValuator bsv = new BlackScholesValuator();
        bsv.setStateOfTheWorld(sotw);
        HashMap<Instrument, Double> prices = new HashMap<Instrument, Double>();
        prices.put(GOOG,488.29);
        prices.put(SP500,1068.76);
        sotw.setPrices(prices);
        HashMap<Option, Double> vols = new HashMap<Option, Double>();
        vols.put(GOOG_OCT_PUT_Option  ,bsv.getImpliedVolatility(GOOG_OCT_PUT_Option  ,(17.3+17.6)/2.0).getPv());
        vols.put(GOOG_DEC_PUT_Option  ,bsv.getImpliedVolatility(GOOG_DEC_PUT_Option  ,(27.5+27.9)/2.0).getPv());

        vols.put(SP500_OCT_PUT_Option ,bsv.getImpliedVolatility(SP500_OCT_PUT_Option ,(25.8+28.2)/2.0).getPv());
        vols.put(SP500_DEC_PUT_Option ,bsv.getImpliedVolatility(SP500_DEC_PUT_Option ,(50.1+52.5)/2.0).getPv());
        sotw.setVolatilities(vols);

        //Set-up a portfolio to iterate over
        List<Option> portfolio = new ArrayList<Option>();
        portfolio.add(GOOG_OCT_PUT_Option  );
        portfolio.add(GOOG_DEC_PUT_Option  );
        portfolio.add(SP500_OCT_PUT_Option );
        portfolio.add(SP500_DEC_PUT_Option );

        int timeSteps=300;
        BinomialTreeValuator btv = new BinomialTreeValuator(sotw,timeSteps);
        TrinomialTreeValuator ttv= new TrinomialTreeValuator(sotw,timeSteps);
        ExplicitFiniteDifferenceValuator efdv= new ExplicitFiniteDifferenceValuator(sotw,timeSteps,timeSteps);
        MonteCarloValuator mcv= new MonteCarloValuator(sotw,timeSteps,100000);
        ImplicitFiniteDifferenceValuator ifdv= new ImplicitFiniteDifferenceValuator(sotw,timeSteps,timeSteps);
        CrankNicholsonFiniteDifferenceValuator cnfdv= new CrankNicholsonFiniteDifferenceValuator(sotw,timeSteps,timeSteps);
        System.out.println("OPTION\tBlack-Scholes\tBinomial\tTrinomial\tEFD\tMonteCarlo\tIFD\tCNFD");
        for (Option option : portfolio){
            System.out.printf("%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",option.getDescription(),
                    bsv.getValuation(option).getPv(),
                    btv.getValuation(option).getPv(),
                    ttv.getValuation(option).getPv(),
                    efdv.getValuation(option).getPv(),
                    mcv.getValuation(option).getPv(),
                    ifdv.getValuation(option).getPv(),
                    cnfdv.getValuation(option).getPv());
        }

    }

    private static void hw2_problem14() throws InvalidInstrumentException, ImmutablesAreImmutableException {
        System.out.println("\n\n===================================================");
        System.out.println("PROBLEM 14");
        Linear testEurUnderlying = new Linear("syntheticUnderlying",null,0.0);
        LocalDate TODAY = new LocalDate(2009,1,1);
        LocalDate EXPIRY = TODAY.plusMonths(6); //T=.5
        Option testEurCallOption = new Option(CALL, AMERICAN, TODAY,100,EXPIRY,testEurUnderlying);
        BaseStateOfTheWorld sotw = new BaseStateOfTheWorld(TODAY,.11);
        HashMap<Instrument, Double> prices = new HashMap<Instrument, Double>();
        HashMap<Option, Double> vols = new HashMap<Option, Double>();
        prices.put(testEurUnderlying,100.0);
        vols.put(testEurCallOption,.3);
        sotw.setPrices(prices);
        sotw.setVolatilities(vols);

        /*
        Calculated value is 10.9571
        For S0=100,  K=100, T=0.5, volatility=0.3, r=0.11, H=110
         */
        int i = 1000000;
        BarrierValuator bmcv = new BarrierValuator(sotw,300,i,110);
        long start = System.currentTimeMillis();
        Valuation mcvValCall = bmcv.getValuation(testEurCallOption);
        long total =  System.currentTimeMillis()-start;
        System.out.println("number of iterations = " + i + " took " + total + "ms");
        System.out.println("Am up and in call is " + mcvValCall.getPv() + " and should be 10.9571");
        System.out.println("bmcv.stddev = " + bmcv.stddev);
        System.out.println("bmcv.stderr = " + bmcv.stderr);
    }
 
    private static void unitTests() throws InvalidInstrumentException, ImmutablesAreImmutableException {
        Linear TEST = new Linear("TEST UNDERLYING",null,0.0);
        LocalDate TEST_DMO = new LocalDate(2010,1,1);
        LocalDate TEST_EXP = new LocalDate(2010,1,1);

        BaseStateOfTheWorld testSotw = new BaseStateOfTheWorld(new LocalDate(2009,1,1),.05);
        Option testCallOption = new Option(CALL, EUROPEAN, TEST_DMO,100,TEST_EXP,TEST);
        Option testPutOption = new Option(PUT, EUROPEAN, TEST_DMO,100,TEST_EXP,TEST);
        Option testCallAmOption = new Option(CALL, AMERICAN, TEST_DMO,100,TEST_EXP,TEST);
        Option testPutAmOption = new Option(PUT, AMERICAN, TEST_DMO,100,TEST_EXP,TEST);
        HashMap<Instrument, Double> pricesTest = new HashMap<Instrument, Double>();
        HashMap<Option, Double> volsTest = new HashMap<Option, Double>();
        pricesTest.put(TEST,100.0);
        volsTest.put(testCallOption,.3);
        volsTest.put(testPutOption,.3);
        volsTest.put(testCallAmOption,.3);
        volsTest.put(testPutAmOption,.3);
        testSotw.setPrices(pricesTest);
        testSotw.setVolatilities(volsTest);

        BlackScholesValuator bsv = new BlackScholesValuator();
        bsv.setStateOfTheWorld(testSotw);
        Valuation testCallVal = bsv.getValuation(testCallOption);
        Valuation testImpliedCallVol = bsv.getImpliedVolatility(testCallOption,14.23124539);
        Valuation testPutVal = bsv.getValuation(testPutOption);
        Valuation testImpliedPutVol = bsv.getImpliedVolatility(testPutOption,9.354187844);
        System.out.printf("calculated\treference\n");
        System.out.printf("%f\t%f\tBlackScholes [%s]\n",testCallVal.getPv(),14.23124539,testCallOption.getDescription());
        System.out.printf("%f\t%f\tBlackScholes [%s]\n",testImpliedCallVol.getPv(),.3,testCallOption.getDescription());
        System.out.printf("%f\t%f\tBlackScholes [%s]\n",testPutVal.getPv(),9.354187844,testPutOption.getDescription());
        System.out.printf("%f\t%f\tBlackScholes [%s]\n",testImpliedPutVol.getPv(),.3,testPutOption.getDescription());

        BinomialTreeValuator btv = new BinomialTreeValuator(testSotw,200);
        btv.setStateOfTheWorld(testSotw);
        testCallVal = btv.getValuation(testCallOption);
        testPutVal = btv.getValuation(testPutOption);
        Valuation testCallAmVal = btv.getValuation(testCallAmOption);
        Valuation testPutAmVal = btv.getValuation(testPutAmOption);
        System.out.printf("calcuated\treference\n");
        System.out.printf("%f\t%f\tBinomialTree [%s]\n",testCallVal.getPv(),14.21653335,testCallOption.getDescription());
        System.out.printf("%f\t%f\tBinomialTree [%s]\n",testPutVal.getPv(),9.339475801,testPutOption.getDescription());
        System.out.printf("%f\t%f\tBinomialTree [%s]\n",testCallAmVal.getPv(),14.21653335,testCallAmOption.getDescription());
        System.out.printf("%f\t%f\tBinomialTree [%s]\n",testPutAmVal.getPv(),9.863161797,testPutAmOption.getDescription());

        TrinomialTreeValuator ttv = new TrinomialTreeValuator(testSotw,2000);
        ttv.setStateOfTheWorld(testSotw);
        testCallVal = ttv.getValuation(testCallOption);
        testPutVal = ttv.getValuation(testPutOption);
        testCallAmVal = ttv.getValuation(testCallAmOption);
        testPutAmVal = ttv.getValuation(testPutAmOption);
        System.out.printf("calcuated\treference\n");
        System.out.printf("%f\t%f\tTrinomialTree [%s]\n",testCallVal.getPv(),14.21653335,testCallOption.getDescription());
        System.out.printf("%f\t%f\tTrinomialTree [%s]\n",testPutVal.getPv(),9.339475801,testPutOption.getDescription());
        System.out.printf("%f\t%f\tTrinomialTree [%s]\n",testCallAmVal.getPv(),14.21653335,testCallAmOption.getDescription());
        System.out.printf("%f\t%f\tTrinomialTree [%s]\n",testPutAmVal.getPv(),9.863161797,testPutAmOption.getDescription());

        ExplicitFiniteDifferenceValuator efdv = new ExplicitFiniteDifferenceValuator(testSotw,2000,2000);
        testCallVal = efdv.getValuation(testCallOption);
        testPutVal = efdv.getValuation(testPutOption);
        testCallAmVal = efdv.getValuation(testCallAmOption);
        testPutAmVal = efdv.getValuation(testPutAmOption);
        System.out.println("\nExplicitFiniteDifferenceValuator");
        System.out.printf("calcuated\treference\n");
        System.out.printf("%f\t%f\tExplicitFiniteDifferenceValuator [%s]\n",testCallVal.getPv(),14.21653335,testCallOption.getDescription());
        System.out.printf("%f\t%f\tExplicitFiniteDifferenceValuator [%s]\n",testPutVal.getPv(),9.339475801,testPutOption.getDescription());
        System.out.printf("%f\t%f\tExplicitFiniteDifferenceValuator [%s]\n",testCallAmVal.getPv(),14.21653335,testCallAmOption.getDescription());
        System.out.printf("%f\t%f\tExplicitFiniteDifferenceValuator [%s]\n",testPutAmVal.getPv(),9.863161797,testPutAmOption.getDescription());

        ImplicitFiniteDifferenceValuator ifdv = new ImplicitFiniteDifferenceValuator(testSotw,2000,2000);
        testCallVal = ifdv.getValuation(testCallOption);
        testPutVal = ifdv.getValuation(testPutOption);
        testCallAmVal = ifdv.getValuation(testCallAmOption);
        testPutAmVal = ifdv.getValuation(testPutAmOption);
        System.out.println("\nImplicitFiniteDifferenceValuator");
        System.out.printf("calcuated\treference\n");
        System.out.printf("%f\t%f\tImplicitFiniteDifference [%s]\n",testCallVal.getPv(),14.21653335,testCallOption.getDescription());
        System.out.printf("%f\t%f\tImplicitFiniteDifference [%s]\n",testPutVal.getPv(),9.339475801,testPutOption.getDescription());
        System.out.printf("%f\t%f\tImplicitFiniteDifference [%s]\n",testCallAmVal.getPv(),14.21653335,testCallAmOption.getDescription());
        System.out.printf("%f\t%f\tImplicitFiniteDifference [%s]\n",testPutAmVal.getPv(),9.863161797,testPutAmOption.getDescription());

        CrankNicholsonFiniteDifferenceValuator cnfdv = new CrankNicholsonFiniteDifferenceValuator(testSotw,2000,2000);
        testCallVal = cnfdv.getValuation(testCallOption);
        testPutVal = cnfdv.getValuation(testPutOption);
        testCallAmVal = cnfdv.getValuation(testCallAmOption);
        testPutAmVal = cnfdv.getValuation(testPutAmOption);
        System.out.println("\nCrankNicholsonFiniteDifferenceValuator");
        System.out.printf("calcuated\treference\n");
        System.out.printf("%f\t%f\tCrankNicholsonFiniteDifferenceValuator [%s]\n",testCallVal.getPv(),14.21653335,testCallOption.getDescription());
        System.out.printf("%f\t%f\tCrankNicholsonFiniteDifferenceValuator [%s]\n",testPutVal.getPv(),9.339475801,testPutOption.getDescription());
        System.out.printf("%f\t%f\tCrankNicholsonFiniteDifferenceValuator [%s]\n",testCallAmVal.getPv(),14.21653335,testCallAmOption.getDescription());
        System.out.printf("%f\t%f\tCrankNicholsonFiniteDifferenceValuator [%s]\n",testPutAmVal.getPv(),9.863161797,testPutAmOption.getDescription());
    }
}
