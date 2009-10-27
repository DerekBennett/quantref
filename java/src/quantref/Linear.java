package quantref;

import org.joda.time.LocalDate;

import java.util.List;
import java.util.ArrayList;

import com.sun.tools.javac.util.Pair;

/**
 * User: derekbennett
 * Date: Oct 22, 2009
 * Time: 11:47:40 AM
 */
public class Linear implements Instrument{
    private String symbol;
    private List<Pair<LocalDate,Double>> dividends;
    private double continuousDividend;

    /**
     * Every linear must have a name. A linear could be a stock, bond, gold lease, swap, etc
     * So, there will be subclasses of Linear in the future.
     * @param symbol
     * @param dividends
     * @param continuousDividend
     */
    public Linear(String symbol, List<Pair<LocalDate, Double>> dividends, Double continuousDividend)
            throws InvalidInstrumentException{
        this.symbol = symbol;
        if (symbol==null){
            throw new InvalidInstrumentException("You must supply a name for every linear instrument");
        }
        if (dividends==null)
            this.dividends=new ArrayList<Pair<LocalDate,Double>>();
        else
            this.dividends = dividends;
        if (continuousDividend==null)
            this.continuousDividend = 0;
        else
            this.continuousDividend = continuousDividend;
    }

    public String getSymbol() {
        return symbol;
    }

    public List<Pair<LocalDate, Double>> getDividends() {
        return dividends;
    }

    public double getContinuousDividend() {
        return continuousDividend;
    }

    public String getDescription() {
        return symbol;
    }
}
