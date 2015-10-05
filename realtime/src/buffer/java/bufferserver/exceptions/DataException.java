package nl.fcdonders.fieldtrip.bufferserver.exceptions;

/**
 * An exception thrown when there is an error/inconsistency within the data.
 *
 * @author wieke
 *
 */
public class DataException extends Exception {

	private static final long serialVersionUID = -4238519569305114389L;

	public DataException(final String string) {
		super(string);
	}

}
