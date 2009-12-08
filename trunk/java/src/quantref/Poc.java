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
    public static Poc parse(String s){
        if (s.equalsIgnoreCase(PUT.toString()))
            return PUT;
        else if (s.equalsIgnoreCase(CALL.toString()))
            return CALL;
        else
            return null;
    }
}
