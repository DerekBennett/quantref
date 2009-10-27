package quantref;

import org.joda.time.LocalDate;

import java.util.Map;

/**
 * User: derekbennett
 * Date: Oct 18, 2009
 * Time: 9:40:48 PM
 * A StateOfTheWorld is all of the data necessary to value a given set of instruments with a given
 * set of Valuators. You may add curves to the environment, but you may not change them, ever. If you want to have
 * the effect of a shifted curve, Created a new implementation of StateOfTheWorld that uses the current one and
 * calculates shifts on-the-fly.
 */
public class BaseStateOfTheWorld implements StateOfTheWorld {
    private LocalDate businessDate;
    private double interestRate;
    private Map<Instrument,Double> prices;
    private Map<Option, Double> volatilities;

    public BaseStateOfTheWorld(LocalDate businessDate, double interestRate) {
        this.businessDate = businessDate;
        this.interestRate = interestRate;
    }

    /**
     * Business date is invariant. If you want to shift the date, create a new environment that delegates to this one
     * @return
     */
    public LocalDate getBusinessDate() {
        return businessDate;
    }

    /**
     * Interest rate is immutable. See business date for instructions.
     * @return
     */
    public double getInterestRate() {
        return interestRate;
    }

    /**
     * Prices are a list of USD prices as of businessDate for the mapped instruments
     * @return
     */
    public Map<Instrument, Double> getPrices() {
        return prices;
    }

    public void setPrices(Map<Instrument, Double> prices) throws ImmutablesAreImmutableException {
        if (this.prices!=null){
            throw new ImmutablesAreImmutableException("You have already set these prices. You may not set them again.");
        }
        this.prices = prices;
    }

    /**
     * Volatilities are a list of volatilities to be used for valuation as of businessDate
     * @return
     */
    public Map<Option, Double> getVolatilities() {
        return volatilities;
    }

    public void setVolatilities(Map<Option, Double> volatilities) throws ImmutablesAreImmutableException {
        if (this.volatilities!=null){
            throw new ImmutablesAreImmutableException("You have already set these volatilities. You may not set them again.");
        }
        this.volatilities = volatilities;
    }
}
