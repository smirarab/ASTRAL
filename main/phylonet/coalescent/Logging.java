package phylonet.coalescent;

import java.util.concurrent.LinkedBlockingQueue;

public class Logging {

	static final String ENDMESSAGE = "&&END&&";
	public static long timer;
	public static final boolean timerOn = false;
	private static LinkedBlockingQueue<String> messageQueue = new LinkedBlockingQueue<String>();
	
	public static void logTimeMessage(String message) {
		if (timerOn) {
			log("TIME TOOK FROM LAST NOTICE " + message 
					+" "+ (double) (System.currentTimeMillis() - timer) / 1000);
			timer = System.currentTimeMillis();
		}
	}
	public static void log(String s) {
		try {
			messageQueue.put(s);
		} catch (InterruptedException e) {
			System.err.println("Message "+s+ " could not be enqueued.");
		}
	}
	
	public static void startLogger() {
		Threading.execute(new Runnable() {

			@Override
			public void run() {
				String message = "";
				while(message!=ENDMESSAGE) {
					System.err.println(message);
					try {
						message = messageQueue.take();
					} catch (InterruptedException e) {
						System.err.println("Message could not be dequeued.");
					}
				}
			}
		});
	}

}
