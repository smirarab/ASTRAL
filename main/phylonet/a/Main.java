import java.io.FileReader;


public class Main {
	public static void main(String[] args) {
		FileReader reader;
		try {
			reader = new FileReader("geneTreesAsInts.txt");
		}
		catch (Exception e) {
			System.out.println("something is wrong :/");
		}
	}
}
