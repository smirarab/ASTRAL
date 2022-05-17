package phylonet.coalescent;

import phylonet.coalescent.Polytree.PTNative;

public class NativeLibraryTester {
	static int[] queue = {-1, 10, 8, 95, 10, 6, 91, 2, 12, 95, 2, 4, 75, 16, 14, 95, 4, 18, 95, 3, 89, 2, 26, 24, 95, 16, 22, 91, 4, 20, 91, 4, 2, 0, 95, 13, 93, 4, -1, 26, 10, 91, 9, 24, 95, 9, 6, 14, 95, 2, 75, 20, 73, 18, 12, 95, 2, 16, 4, 95, 7, 8, 75, 22, 38, 95, 3, 75, 73, -1, 36, 10, 91, 6, 8, 95, 11, 16, 95, 5, 6, 89, 2, 20, 14, 95, 3, 12, 79, 0, 2, 22, 4, 91, 2, 18, 95, 2, 75, 79, 89, 2, -1, 8, 4, 91, 2, 42, 91, 2, 6, 95, 3, 12, 89, 2, 22, 32, 20, 18, 38, 95, 3, 75, 75, 73, -1, 24, 10, 95, 4, 8, 91, 3, 4, 91, 3, 16, 91, 3, 26, 93, 3, 18, 6, 95, 8, 20, 22, 2, 14, 91, 3, 0, 12, 91, 3, 95, 3, 95, 2, 75, 89, 2, -1, 54, 14, 4, 91, 3, 89, 3, 22, 12, 2, 75, 72, 75, 20, 75, 0, 91, 3, 89, 3, -1, 52, 18, 22, 95, 4, 95, 2, 12, 73, 48, 6, 91, 3, 20, 91, 3, 2, 91, 3, 0, 95, 3, 14, 73, -1, 28, 72, 95, 3, 2, 75, 4, 75, 20, 73, 36, 14, 91, 2, 12, 91, 2, 16, 22, 91, 3, 91, 3, 0, 89, 3, -1, 36, 8, 91, 3, 22, 91, 3, 4, 10, 91, 3, 16, 91, 3, 89, 3, 2, 18, 20, 6, 95, 3, 14, 75, 12, 95, 3, 0, 95, 3, 75, 73, -1, 28, 24, 91, 3, 26, 91, 3, 18, 91, 2, 22, 91, 2, 20, 89, 3, 6, 4, 95, 3, 16, 38, 14, 12, 95, 5, 91, 2, 75, 73, -1, 28, 12, 91, 2, 4, 95, 2, 34, 73, 6, 40, 73, -1, 68, 26, 73, 52, 48, 91, 2, 20, 73, 6, 18, 14, 95, 2, 12, 75, 50, 75, 73, -1, 42, 16, 75, 8, 73, 12, 20, 91, 2, 14, 89, 2, 2, 0, 60, 91, 2, 89, 2, -1, 42, 4, 75, 8, 73, 64, 14, 73, 20, 16, 75, 22, 12, 75, 66, 75, 73, -1, 22, 20, 74, 75, 73, -1, 46, 86, 75, 2, 89, 2, -1, 52, 78, 14, 75, 73, 12, 82, 73, -1, 84, 4, 91, 2, 20, 91, 2, 2, 89, 2, 36, 12, 75, 14, 73, -1, 18, 2, 94, 20, 75, 6, 73, 90, 91, 2, 89, 2, -1, 52, 22, 75, 18, 73, 16, 92, 74, 75, 73, -1, 34, 4, 75, 30, 73, -1, 70, 18, 73, 44, 20, 73, 88, 50, 73, -1, 54, 58, 73, 6, 78, 4, 73, 62, 73, -1, 56, 22, 75, 16, 66, 75, 73, -1, 72, 20, 95, 2, 76, 73, -1, 100, 12, 73, -1, 80, 94, 73, -1, -1, 56, 6, 73, -1, 92, 16, 75, 0, 2, 94, 75, 75, 73, -1, 18, 16, 75, 14, 73, 96, 6, 73, -1, 86, 12, 75, 98, 75, 22, 73, 48, 52, 38, 75, 73};
	static long[][] bit = {{16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}, {16383}};
	static long[][] a = {{13108}, {520}, {12348}, {12848}, {1024}, {2243}, {520}, {13108}, {1024}, {14839}, {2243}, {2243}, {192}, {15807}, {13756}, {2051}, {512}, {15679}, {576}, {512}, {704}, {13480}, {1544}, {1096}, {520}, {576}, {14647}, {512}, {520}, {704}, {14647}, {14647}};
	static long[][] b = {{8}, {12596}, {3267}, {3267}, {2243}, {512}, {2243}, {2243}, {14839}, {512}, {12596}, {1292}, {14140}, {512}, {2051}, {1216}, {15743}, {640}, {15679}, {192}, {2051}, {2327}, {14775}, {512}, {1088}, {14775}, {1608}, {14647}, {14647}, {14647}, {1160}, {192}};
	static long[][] c = {{3267}, {3267}, {768}, {268}, {13116}, {13628}, {13620}, {1032}, {520}, {1032}, {1544}, {12848}, {2051}, {64}, {576}, {13116}, {128}, {64}, {128}, {15679}, {13628}, {576}, {64}, {14775}, {14775}, {1032}, {128}, {1224}, {1216}, {1032}, {576}, {1544}};
	static long[] result = new long[32];
	public static void main(String[] args) {
		
		Factory.instance = new FactoryAstralMP();
		
		Logging.initalize(Factory.instance.newLogger());
		Threading.startThreading(2);
		Logging.startLogger();
		// TODO Auto-generated method stub
		boolean useNativeMethod = true;
		try {
			System.loadLibrary("Astral");
			Logging.log("Native AVX library found.");
			useNativeMethod = true;
		}
		catch (Throwable e) {
			useNativeMethod = false;
			Logging.log("Fail to load native library "+System.mapLibraryName("Astral")+"! Is library path set correctly using -Djava.library.path=lib/?");
			Logging.log("\n\n" + e);
		}
		
		if (useNativeMethod) {
			try {
				PTNative.cppInit(14, 51, queue, bit);
				PTNative.cppBatchCompute(result, a, b, c);
				Logging.log("Native AVX library functions correctly as expected.");
			}
			catch (Throwable e) {
				Logging.log("Native AVX library does not function correctly!");
				Logging.log("Please run compile_native_c_code.sh to compile the library for your machine.");
				Logging.log("Or you can compile the library using following command:");
				Logging.log("g++ -std=c++11 -I\"PATH_TO_FOLDER_CONTAINING_JNI_DOT_H\" -I\"PATH_TO_FOLDER_CONTAINING_JNI_DOT_H/OS\" -march=native -Ofast -fPIC -o lib/libAstral.so -shared main/phylonet_coalescent_Polytree_PTNative.cpp");
				Logging.log("You can also set -D\"java.library.path=lib/no_avx2\" to check whether such version of library functions correctly.");
				Logging.log("\n\n" + e);
			}
		}
		
		Logging.endLogger();
		Threading.shutdown();
	}

}
