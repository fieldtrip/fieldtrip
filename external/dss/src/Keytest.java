import java.io.IOException;

public class Keytest {
    public static boolean test(char key) {
	try {
	    if (System.in.available()==0) return false;
	    return System.in.read()==(int)key;
	} catch (IOException e) {
	    return false;
	}
    }
}
