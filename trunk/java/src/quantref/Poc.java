package quantref;

/**
 * User: derekbennett
 * Date: Oct 24, 2009
 * Time: 12:19:37 PM
 */
public enum Poc {
    PUT("put"),CALL("call");
    private String s;

    Poc(String s) {
        this.s = s;
    }

    @Override
    public String toString() {
        return s;
    }
}
