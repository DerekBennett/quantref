package quantref;

/**
 * User: derekbennett
 * Date: Oct 18, 2009
 * Time: 9:51:50 PM
 */
public class Valuation<T> {
    private T instrument;
    private double pv;

    public Valuation(T instrument, double pv) {
        this.instrument = instrument;
        this.pv = pv;
    }
    public T getInstrument(){
        return instrument;
    }
    public double getPv(){
        return pv;
    }
}
