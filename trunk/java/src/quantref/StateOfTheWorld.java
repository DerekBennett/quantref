package quantref;

import org.joda.time.LocalDate;

import java.util.Map;

/**
 * User: derekbennett
 * Date: Oct 24, 2009
 * Time: 12:15:09 PM
 */
public interface StateOfTheWorld {
    public LocalDate getBusinessDate();
    public double getInterestRate();
    public Map<Instrument, Double> getPrices();
    public Map<Option, Double> getVolatilities();
}
