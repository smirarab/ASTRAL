package phylonet.coalescent;

public class Logging {

	private static long timer;
	private static final boolean timerOn = false;
	private static LoggerInterface logger = null;
	
	public static void logTimeMessage(String message) {
		if (timerOn) {
			logger.log("TIME TOOK FROM LAST NOTICE " + message 
					+" "+ (double) (System.currentTimeMillis() - timer) / 1000);
			timer = System.currentTimeMillis();
		}
	}
	public static void log(String s) {
		logger.log(s);
	}
	
	public static void startLogger() {
		logger.startLogger();
		if (Logging.timerOn) {
			Logging.log("Timer starts here");
			Logging.timer = System.currentTimeMillis();
		}
	}
	
	public static void endLogger() {
		logger.endLogger();
	}
	public static void initalize(LoggerInterface theLogger) {
		logger = theLogger;
		
	}

}
