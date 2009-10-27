package quantref;

import org.joda.time.LocalDate;


/**
 * User: derekbennett
 * Date: Oct 18, 2009
 * Time: 9:49:22 PM
 * Immutable representing all the information you need to know about an option except to value it given an
 * environment and a valuator
 */
public class Option implements Instrument {
    private Poc poc;
    private OptionType type;
    private LocalDate dmo;
    private double strike;
    private LocalDate expiry;
    private Linear underlying;

    public Option(Poc poc, OptionType type, LocalDate dmo, double strike, LocalDate expiry, Linear underlying)
            throws InvalidInstrumentException {
        if (underlying==null) {
            throw new InvalidInstrumentException("Underlying was null when constructing option[" +
                    dmo + " " + poc + " @" + strike + "]");
        }
        this.underlying = underlying;
        if (poc==null) {
            throw new InvalidInstrumentException("Put Or Call indicator was null when constructing option for[" +
                    underlying.getDescription() + "]");
        }
        this.poc = poc;
        if (type==null) {
            throw new InvalidInstrumentException("Option type indicator was null when constructing option for[" +
                    underlying.getDescription() + "]");
        }
        this.type = type;
        if (dmo==null) {
            throw new InvalidInstrumentException("DMO indicator was null when constructing option for[" +
                    underlying.getDescription() + " " + poc + " @" + strike + "]");
        }
        this.dmo = dmo;
        if (expiry==null) {
            throw new InvalidInstrumentException("Expiry was null when constructing option[" + getDescription() + "]");
        }
        this.expiry = expiry;
        if (strike<=0.0) {
            throw new InvalidInstrumentException("Strike was zero or less when constructing option[" + getDescription() + "]");
        }
        this.strike = strike;
    }

    public Poc getPoc() {
        return poc;
    }

    public OptionType getType() {
        return type;
    }

    public LocalDate getDmo() {
        return dmo;
    }

    public double getStrike() {
        return strike;
    }

    public LocalDate getExpiry() {
        return expiry;
    }

    public Linear getUnderlying() {
        return underlying;
    }

    public String getDescription() {
        return underlying.getDescription() + " " + dmo + " " + type + " " + poc + " @" + strike;
    }

    public boolean isPut() {
        return getPoc().equals(Poc.PUT);
    }
    public boolean isCall() {
        return getPoc().equals(Poc.CALL);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Option option = (Option) o;

        if (Double.compare(option.strike, strike) != 0) return false;
        if (!dmo.equals(option.dmo)) return false;
        if (!expiry.equals(option.expiry)) return false;
        if (poc != option.poc) return false;
        if (type != option.type) return false;
        if (!underlying.equals(option.underlying)) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = poc.hashCode();
        result = 31 * result + type.hashCode();
        result = 31 * result + dmo.hashCode();
        temp = strike != +0.0d ? Double.doubleToLongBits(strike) : 0L;
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        result = 31 * result + expiry.hashCode();
        result = 31 * result + underlying.hashCode();
        return result;
    }
}
