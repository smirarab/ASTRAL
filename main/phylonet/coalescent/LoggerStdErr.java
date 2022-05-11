package phylonet.coalescent;

public class LoggerStdErr implements LoggerInterface  {
	
	public void log(String s) {
			System.err.println(s);
	}
	
	@Override
	public void startLogger() {
	}

	@Override
	public void endLogger() {
		
	}

}
