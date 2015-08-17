package nl.fcdonders.fieldtrip.bufferserver.exceptions;

/**
 * An exception thrown when the client does something it shouldn't.
 *
 * @author wieke
 *
 */
public class ClientException extends Exception {

	private static final long serialVersionUID = -7945125817344034976L;

	public ClientException(final String string) {
		super(string);
	}
}
