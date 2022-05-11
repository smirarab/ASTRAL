package phylonet.coalescent;

import java.util.concurrent.LinkedBlockingQueue;

public class LoggerThreaded implements LoggerInterface {

	static final String ENDMESSAGE = "&&END&&";

	private static LinkedBlockingQueue<String> messageQueue = new LinkedBlockingQueue<String>();

	@Override
	public void log(String s) {
		try {
			messageQueue.put(s);
		} catch (InterruptedException e) {
			System.err.println("Message "+s+ " could not be enqueued.");
		}
	}

	@Override
	public void endLogger() {
		this.log(ENDMESSAGE);
	}
	@Override
	public void startLogger() {
		Threading.execute(new Runnable() {

			@Override
			public void run() {
				String message = "";
				while(true) {
					try {
						message = messageQueue.take();
						if (message==ENDMESSAGE)
							break;
						System.err.println(message);
					} catch (InterruptedException e) {
						System.err.println("Message could not be dequeued.");
					}
				}
			}
		});
	}
}
