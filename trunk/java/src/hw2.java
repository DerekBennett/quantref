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
/*
        File data1 = new File("/Users/derekbennett/Documents/stevens/FE621/hw1/data1.csv");
        System.out.println("data1.canRead() = " + data1.canRead());
*/

    }
}
