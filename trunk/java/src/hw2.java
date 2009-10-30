import quantref.*;
import static quantref.OptionType.*;
import static quantref.Poc.*;

import java.util.HashMap;

import org.joda.time.LocalDate;

/**
 * User: derekbennett
 * Date: Oct 18, 2009
 * Time: 2:36:26 PM
 */
public class hw2 {
    public static void main(String[] args) throws InvalidInstrumentException, ImmutablesAreImmutableException {
        Linear GOOG = new Linear("GOOG",null,0.0);
        Linear SP500 = new Linear("SP500",null,0.0);
        Linear TEST = new Linear("TEST UNDERLYING",null,0.0);

        LocalDate OCT09_DMO = new LocalDate(2009,10,1);
        LocalDate DEC09_DMO = new LocalDate(2009,12,1);
        LocalDate TEST_DMO = new LocalDate(2010,1,1);

        LocalDate OCT09_EXP = new LocalDate(2009,10,16);
        LocalDate DEC09_EXP = new LocalDate(2009,12,18);
        LocalDate TEST_EXP = new LocalDate(2010,1,1);

        BaseStateOfTheWorld sotw1 = new BaseStateOfTheWorld(new LocalDate(2009,9,16),.04);
        HashMap<Instrument, Double> prices1 = new HashMap<Instrument, Double>();
        prices1.put(GOOG,488.29);
        prices1.put(SP500,1068.76);
        sotw1.setPrices(prices1);

        BaseStateOfTheWorld sotw = new BaseStateOfTheWorld(new LocalDate(2009,9,17),.04);
        HashMap<Instrument, Double> prices2 = new HashMap<Instrument, Double>();
        prices2.put(GOOG,491.72);
        prices2.put(SP500,1065.49);
        sotw.setPrices(prices2);

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


        hw2_problem1();
        hw2_problem2();
        hw2_problem3();
        hw2_problem4();
/*
        File data1 = new File("/Users/derekbennett/Documents/stevens/FE621/hw1/data1.csv");
        System.out.println("data1.canRead() = " + data1.canRead());
*/

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
        Valuation btVal = btv.getValuation(synthPutAmOption);

        TrinomialTreeValuator ttv = new TrinomialTreeValuator(sotw,200);
        Valuation ttVal = ttv.getValuation(synthPutAmOption);

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
        System.out.println("Black-Scholes price  = " + bsVal.getPv());
        System.out.println("Binomial Tree price  = " + btVal.getPv());
        System.out.println("Trinomial Tree price = " + ttVal.getPv());
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
        System.out.println("PROBLEM 1");
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
    
    private static void hw2_problem4() {
        System.out.println("\n\n===================================================");
        System.out.println("PROBLEM 1");
        System.out.println("Implement an explicit Þnite di?erence method for the same purpose as the \n" +
                "part (1).\n");
        System.out.println("ANSWER\n");
    }
}
