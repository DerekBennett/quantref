import au.com.bytecode.opencsv.CSVReader;
import au.com.bytecode.opencsv.CSVWriter;

import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.HashMap;

import quantref.*;
import static quantref.OptionType.EUROPEAN;
import org.joda.time.LocalDate;
import com.sun.tools.javac.util.Pair;

/**
 * User: derekbennett
 * Date: Dec 6, 2009
 * Time: 10:02:44 AM
 */
public class hw4 {
    private static final String BASEDIR = "/Users/derekbennett/Documents/stevens/FE621-Computational Methods in Finance/hw4/";

    public static void main(String[] args) throws IOException, InvalidInstrumentException, ImmutablesAreImmutableException {
        problem1();
        problem2();
        problem3();
        problem4();

        bonus1();
        bonus2();
        bonus3();
    }

    private static void problem1() {
    }
    private static void problem2() throws IOException, InvalidInstrumentException, ImmutablesAreImmutableException {
        // set-up the world
        final LocalDate VALUE_DATE = new LocalDate(2009,11,23);
        final LocalDate JAN10 = new LocalDate(2010,1,1);
        final LocalDate JAN10_expiry = new LocalDate(2009,12,18);
        final Linear GOOG = new Linear("GOOG",new ArrayList<Pair<LocalDate, Double>>(),0.0 );
        // Federal funds (effective) rate as of 11/23/2009
        final double r=.0012; // http://www.federalreserve.gov/releases/H15/data/Daily/H15_FF_O.txt
        final BaseStateOfTheWorld bsotw = new BaseStateOfTheWorld(VALUE_DATE,r);
        final HashMap<Instrument, Double> prices = new HashMap<Instrument, Double>();
        prices.put(GOOG,582.35);
        bsotw.setPrices(prices);
        BlackScholesValuator bsv = new BlackScholesValuator();

        // read CSV
        FileReader fileReader = null;
        fileReader = new FileReader(BASEDIR + "GOOGoptiondataNov23price=582.35.csv");
        CSVReader csv = new CSVReader(fileReader);
        List<String[]> allRows = csv.readAll();
        
        // loop through rows and calculate implied volatilities
        ArrayList<Option> options = new ArrayList<Option>();
        ArrayList<Valuation> impVols = new ArrayList<Valuation>();
        ArrayList<Valuation> marks = new ArrayList<Valuation>();
        for (String[] row : allRows){
            if (row[0].equals("Poc")) continue; // skip header
            try {
                Option o = new Option(Poc.parse(row[0]), EUROPEAN,JAN10,Double.parseDouble(row[1]),JAN10_expiry,GOOG);
                options.add(o);
                bsv.setStateOfTheWorld(bsotw);
                // Use average of BID and ASK as price from which to imply vols for each option
                Valuation m = new Valuation<Option>(o,(Double.parseDouble(row[5]) + Double.parseDouble(row[6])) / 2.0);
                Valuation v = bsv.getImpliedVolatility(o, m.getPv());
                if (v.getPv()>.99 || v.getPv()<.01){
                    // If vols still blow-up, exclude them
                    System.out.printf("row with poc,strike=[%s,%s] had a bad vols and was skipped\n",row[0],row[1]);
                    continue;
                }
                impVols.add(v);
                marks.add(m);
            }
            catch (NumberFormatException nfe){
                System.out.printf("row with poc,strike=[%s,%s] had a bad bid or ask and was skipped\n",row[0],row[1]);
                continue;
            }
        }

        // Output results for graphing
        CSVWriter writer = new CSVWriter(new FileWriter(BASEDIR+"hw4-problem2-implied-vols.csv"));
        for(Valuation v : impVols){
            String [] line = new String[2];
            line[0] = String.valueOf(((Option)v.getInstrument()).getStrike());
            line[1] = String.valueOf(v.getPv());
            writer.writeNext(line);
        }
        writer.flush();
        writer.close();

    }
    private static void problem3() {
    }
    private static void problem4() {
    }
    private static void bonus1() {
    }
    private static void bonus2() {
    }
    private static void bonus3() {
    }
}
