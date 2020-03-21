package phylonet.coalescent;

public class Logging {

	public static long timer;
	public static final boolean timerOn = true;
	public static void logTimeMessage(String message) {
		if (timerOn) {
			System.err.println("TIME TOOK FROM LAST NOTICE " + message 
					+" "+ (double) (System.currentTimeMillis() - timer) / 1000);
			timer = System.currentTimeMillis();
		}
	}

}
