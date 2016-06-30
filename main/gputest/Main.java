import java.io.FileReader;
import java.util.ArrayList;
import java.util.Scanner;


public class Main {
	public static final int SPECIES_LENGTH = 200;
	public static void main(String[] args) {
		ArrayList<Integer> gtaiArrayList= new ArrayList<Integer>();
		ArrayList<int[]> allArrayList = new ArrayList<int[]>();
		ArrayList<int[]> tripartitions = new ArrayList<int[]>();
		ArrayList<Integer> weight = new ArrayList<Integer>();
		FileReader reader;
		try {
			reader = new FileReader("geneTreesAsInts.txt");
			Scanner in = new Scanner(reader);
			while(in.hasNextInt()){
				gtaiArrayList.add(in.nextInt());
			}
			in.close();
			reader.close();
			reader = new FileReader("Tripartitions.txt");
			int read;
			String integerString = "";
			int[] allArray = new int[SPECIES_LENGTH];
			while((read = reader.read()) != '&') {
				if(read == '}') {
					allArray[Integer.parseInt(integerString)] = 1;
					integerString = "";
					allArrayList.add(allArray.clone());
					allArray = new int[SPECIES_LENGTH];
				}
				if(read <= '9' && read >= '0') {
					integerString += (char)read;
				}
				if(read == ',') {
					allArray[Integer.parseInt(integerString)] = 1;
					integerString = "";
				}
			}
			int tripCounter = 1;
			in = new Scanner(reader);
			while((read = reader.read()) != -1) {
				if(read == '}') {
					allArray[Integer.parseInt(integerString)] = 1;
					integerString = "";
					tripartitions.add(allArray.clone());
					allArray = new int[SPECIES_LENGTH];
					tripCounter++;
					if(tripCounter == 4) {
						tripCounter = 1;
						reader.read();
						reader.read();
						while((read = reader.read()) != -1 && read != '\r') {
							integerString += (char)read;
						}
						weight.add(Integer.parseInt(integerString));
						integerString = "";
					}
				}
				if(read <= '9' && read >= '0') {
					integerString += (char)read;
				}
				if(read == ',') {
					allArray[Integer.parseInt(integerString)] = 1;
					integerString = "";
				}
			}

		}
		catch (Exception e) {
			e.printStackTrace();
			System.out.println("something is wrong :/");
		}
		int[] geneTreesAsInts = new int[gtaiArrayList.size()];
		for(int i = 0; i < geneTreesAsInts.length; i++) {
			geneTreesAsInts[i] = gtaiArrayList.get(i);
		}
		int[] allArray1D = new int[allArrayList.size() * Main.SPECIES_LENGTH];
		int[] tripartitions1D = new int[tripartitions.size() * Main.SPECIES_LENGTH];
		for(int i = 0; i < allArrayList.size(); i++) {
			for(int j = 0; j < SPECIES_LENGTH; j++) {
				allArray1D[i * SPECIES_LENGTH + j] = allArrayList.get(i)[j];
			}
		}
		for(int i = 0; i < tripartitions.size(); i++) {
			for(int j = 0; j < SPECIES_LENGTH; j++) {
				tripartitions1D[i * SPECIES_LENGTH + j] = tripartitions.get(i)[j];

			}
		}
		GPUCall call = new GPUCall(geneTreesAsInts, allArray1D, tripartitions1D);

		long timeA;
		for(int i = 12; i <= 13; i++) {
			call.WORK_GROUP_SIZE = 1L << i;		
			timeA = System.nanoTime();
			call.update();
			System.out.println(System.nanoTime()-timeA + " " + i);
		}	
		
	}
}
