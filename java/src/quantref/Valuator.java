package quantref;

/**
 * User: derekbennett
 * Date: Oct 18, 2009
 * Time: 10:02:38 PM
 */
public interface Valuator<T extends Instrument> {
    public void setStateOfTheWorld(final StateOfTheWorld sotw);
    public Valuation getValuation(final T instrument);
    public Valuation getImpliedVolatility(final T instrument, final double price);
}
